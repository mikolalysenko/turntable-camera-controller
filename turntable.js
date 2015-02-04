'use strict'

module.exports = createTurntableController

var filterVector = require('filtered-vector')
var lookAt       = require('gl-mat4/lookAt')
var cross        = require('gl-vec3/cross')
var normalize3   = require('gl-vec3/normalize')
var dot3         = require('gl-vec3/dot')

function len3(x, y, z) {
  return Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2))
}

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

function TurntableController(zoomMin, zoomMax, center, up, right, radius, theta, phi) {
  this.center = filterVector(center)
  this.up     = filterVector(up)
  this.right  = filterVector(right)
  this.radius = filterVector([radius])
  this.angle  = filterVector([theta, phi])
  this.angle.bounds = [[-Infinity,-Math.PI/2], [Infinity,Math.PI/2]]
  this.setZoomBounds(zoomMin, zoomMax)

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
  this._lastTick      = 0

  this._recalcMatrix(0)
}

var proto = TurntableController.prototype

proto.dirty = function() {
  return this._isDirty
}

proto.get = function(result) {
  this._isDirty = false
  if(!result) {
    return this.computedMatrix
  }
  var m = this.computedMatrix
  for(var i=0; i<16; ++i) {
    result[i] = m[i]
  }
  return result
}

proto.setZoomBounds = function(minDist, maxDist) {
  minDist = minDist || 0.0
  maxDist = maxDist || Infinity
  this._isDirty = true
}

proto._recalcMatrix = function(t) {
  //Recompute relevant curves
  this.center.curve(t)
  this.up.curve(t)
  this.right.curve(t)
  this.radius.curve(t)
  this.angle.curve(t)

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

  //Compute toward vector
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
    this.right.lastT(),
    this.radius.lastT(),
    this.angle.lastT())
  if(t1 < t0 && 
    this.center.stable() &&
    this.up.stable() &&
    this.right.stable() &&
    this.radius.stable() &&
    this.angle.stable()) {
    return
  }
  this._recalcMatrix(t)
  this._lastTick = t
  this._isDirty = true
}

proto.rotate = function(t, dtheta, dphi) {
  var radius = Math.exp(this.computedRadius[0])
  this.angle.move(t, dtheta / radius, dphi / radius)
}

proto.pan = function(t, dx, dy) {
  var mat = this.computedMatrix
  var vx = mat[0] * dx + mat[4] * dy
  var vy = mat[1] * dx + mat[5] * dy
  var vy = mat[2] * dx + mat[6] * dy
  this.center.move(t, vx, vy, vz)
}

//Recenters the 
proto.tare = function(t) {
  //Recompute state for new t value
  var prevT = this._lastTick
  this._recalcMatrix(t)

  var mat = this.computedMatrix

  var ux = mat[0]
  var uy = mat[1]
  var uz = mat[2]
  var ul = len3(ux, uy, uz)
  ux /= ul
  uy /= ul
  uz /= ul

  var rx = mat[4]
  var ry = mat[5]
  var rz = mat[6]
  var rl = len3(rx, ry, rz)
  rx /= rl
  ry /= rl
  rz /= rl

  this.up.jump(t, ux, uy, uz)
  this.right.jump(t, rx, ry, rz)
  this.angle.jump(t, 0, 0)

  //Reset state of coordinates
  this._recalcMatrix(prevT)
}

proto.idle = function(t) {
  this.center.idle(t)
  this.up.idle(t)
  this.right.idle(t)
  this.radius.idle(t)
  this.angle.idle(t)
}

proto.flush = function(t) {
  this.center.flush(t)
  this.up.flush(t)
  this.right.flush(t)
  this.radius.flush(t)
  this.angle.flush(t)
}

proto.lookAt = function(t, center, up, eye) {
  center = center || this.computedCenter
  up     = up     || this.computedUp
  eye    = eye    || this.computedEye

  var ux = up[0]
  var uy = up[1]
  var uz = up[2]
  var ul = len3(ux, uy, uz)
  if(ul < 1e-6) {
    return
  }
  ux /= ul
  uy /= ul
  uz /= ul

  var tx = eye[0] - center[0]
  var ty = eye[1] - center[1]
  var tz = eye[2] - center[2]
  var tl = len3(tx, ty, tz)
  if(tl < 1e-6) {
    return
  }
  tx /= tl
  ty /= tl
  tz /= tl

  var right = this.computedRight
  var rx = right[0]
  var ry = right[1]
  var rz = right[2]
  var ru = ux*rx + uy*ry + uz*rz
  rx -= ru * ux
  ry -= ru * uy
  rz -= ru * uz
  var rl = len3(rx, ry, rz)

  if(rl > 1e-6) {
    rx = uy * tz - uz * ty
    ry = uz * tx - ux * tz
    rz = ux * ty - uy * tx
    rl = len3(rx, ry, rz)
    if(rl < 1e-6) {
      return
    }
  }
  rx /= rl
  ry /= rl
  rz /= rl

  this.up.set(t, ux, uy, uz)
  this.right.set(t, rx, ry, rz)
  this.center.set(t, center[0], center[1], center[2])

  var tu = ux*tx + uy*ty + uz*tz
  var ru = rx*tx + ry*ty + rz*tz

  var phi = Math.acos(tu)
  var theta = Math.acos(ru)

  this.radius.set(t, Math.log(tl))
  this.angle.set(t, theta, phi) //TODO: Prevent coordinate wrap around here
}

function createTurntableController(options) {
  options = options || {}
  var center = options.center || [0,0,0]
  var up     = options.up     || [0,1,0]
  var right  = options.right  || findOrthoPair(up)
  var radius = options.radius || 1.0
  var theta  = options.theta  || 0.0
  var phi    = options.phi    || 0.0

  center = [].slice.call(center, 0, 3)

  up = [].slice.call(up, 0, 3)
  normalize3(up, up)

  right = [].slice.call(right, 0, 3)
  normalize3(right, right)

  if('eye' in options) {
    var eye = options.eye
    var toward = [
      eye[0]-center[0],
      eye[1]-center[1],
      eye[2]-center[2]
    ]
    cross(right, toward, up)
    if(len3(right[0], right[1], right[2]) < 1e-6) {
      right = findOrthoPair(up)
    } else {
      normalize3(right, right)
    }

    radius = len3(toward[0], toward[1], toward[2])

    var ut = dot3(up, toward) / radius
    var rt = dot3(right, toward) / radius
    phi   = Math.acos(ut)
    theta = Math.acos(rt)
  }

  //Use logarithmic coordinates for radius
  radius = Math.log(radius)

  //Return the controller
  return new TurntableController(
    options.zoomMin,
    options.zoomMax,
    center,
    up,
    right,
    radius,
    theta,
    phi)
}