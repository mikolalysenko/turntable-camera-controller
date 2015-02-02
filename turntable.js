'use strict'

var filterVector = require('filtered-vector')
var lookAt = require('gl-mat4/look-at')
var cross  = require('gl-vec3/cross')
var normalize3 = require('gl-vec3/normalize')

function findOrthoPair(v) {
  var vx = Math.abs(v[0])
  var vy = Math.abs(v[1])
  var vz = Math.abs(v[2])

  var u = [0,0,0]
  if(vx > Math.max(vy, vz)) {
    u[2] = 1
  } else if(vy > Math.max(vx, vz)) {
    u[0] = 1
  } else {
    u[1] = 1
  }

  var vv = 0
  var uv = 0
  for(var i=0; i<3; ++i ) {
    vv += v[i] * v[i]
    uv += u[i] * v[i]
  }
  for(var i=0; i<3; ++i) {
    u[i] -= (uv / vv) *  v[i]
  }
  normalize3(u, u)
  return u
}

function TurntableController(center, up, radius, theta, phi) {
  this.center = filterVector(center)
  this.up     = filterVector(up)
  this.right  = filterVector(findOrthoPair(up))
  this.radius = filterVector([radius])
  this.angle  = filterVector([theta, phi])

  this.computedCenter = this.center.curve(0)
  this.computedUp     = this.up.curve(0)
  this.computedRight  = this.right.curve(0)
  this.computedRadius = this.radius.curve(0)
  this.computedAngle  = this.angle.curve(0)
  this.computedToward = [0,0,0]
  this.computedEye    = [0,0,0]
  this.computedMatrix = new Array(16)

  for(var i=0; i<16; ++i) {
    this.computedMatrix[i] = 0.0
  }

  this._isDirty       = true
  this._lastTick      = -Infinity

  this.recalcMatrix()
}

var proto = TurntableController.prototype

proto.dirty = function() {
  return this._isDirty
}

proto.get = function(result) {
  var m = this.computedMatrix
  for(var i=0; i<16; ++i) {
    result[i] = m[i]
  }
  this._isDirty = false
}

proto._recalcMatrix = function() {

  //Compute frame for camera matrix
  var up     = this.computedUp
  var right  = this.computedRight
  var uu = 0.0
  var ur = 0.0
  for(var i=0; i<3; ++i) {
    ur += up[i] * right[i]
    uu += up[i] * up[i]
  }
  var ul = Math.sqrt(uu)
  var rr = 0.0
  for(var i=0; i<3; ++i) {
    right[i] -= up[i] * ur / uu
    rr       += right[i] * right[i]
    up[i]    /= ul
  }
  var rl = Math.sqrt(rr)
  for(var i=0; i<3; ++i) {
    right[i] /= rl
  }
  var toward = this.computedToward
  cross(toward, up, right)
  normalize3(toward, toward)

  //Compute angular parameters
  var radius = Math.exp(this.computedRadius[0])
  var theta  = this.computedAngle[0] % (2.0 * Math.PI)
  var phi    = this.computedAngle[1]

  var ctheta = Math.cos(theta)
  var stheta = Math.sin(theta)
  var cphi   = Math.cos(phi)
  var sphi   = Math.sin(phi)

  var center = this.computedCenter

  var wx = radius * ctheta * cphi 
  var wy = radius * stheta * cphi
  var wz = radius * sphi

  var eye = this.computedEye
  for(var i=0; i<3; ++i) {
    eye[i] = center[i] + wx * right[i] + wy * toward[i] + wz * up[i]
  }

  lookAt(
    this.computedMatrix,
    eye,
    this.computedCenter,
    this.computedUp)
}

proto.tick = function(t) {
  var t0 = Math.min(this._lastTick, t)
  var t1 = Math.max(
    this.center.lastT(),
    this.up.lastT(),
    this.radius.lastT(),
    this.angle.lastT())
  if(t1 < t0 && 
    this.center.stable() &&
    this.up.stable() &&
    this.radius.stable() &&
    this.angle.stable()) {
    return
  }
  this.center.curve(t)
  this.up.curve(t)
  this.radius.curve(t)
  this.angle.curve(t)
  this._recalcMatrix()
  this._lastTick = t
  this._isDirty = true
}

proto.flush = function(t) {
  this.center.flush(t)
  this.up.flush(t)
  this.radius.flush(t)
  this.angle.flush(t)
}

proto.zoom = function(t, dr) {
  this.radius.move(t, dr)
}

proto.translate = function(t, dx, dy, dz) {
  this.center.move(t, dx, dy, dz)
}

proto.rotate = function(t, dtheta, dphi) {
  this.angle.move(t, dtheta, dphi)

  //hack:  clamp phi to range [-pi/2, pi/2]
  var state = this.angle._state
  var n     = state.length - 2
  state[n] = Math.max(Math.min(state[n], Math.PI/2), Math.PI/2)
}

proto.pan = function(t, dx, dy) {
  var mat = this.computedMatrix
  var vx = mat[0] * dx + mat[4] * dy
  var vy = mat[1] * dx + mat[5] * dy
  var vy = mat[2] * dx + mat[6] * dy
  this.center.move(t, vx, vy, vz)
}

function createTurntableController(options) {
  options = options || {}
  return new TurntableController(
    options.center || [0,0,0],
    options.up     || [0,1,0],
    options.right  || [1,0,0],
    options.radius || 1.0,
    options.theta  || 0.0,
    options.phi    || 0.0)
}