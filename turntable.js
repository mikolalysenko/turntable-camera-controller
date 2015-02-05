'use strict'

module.exports = createTurntableController

var filterVector = require('filtered-vector')
var invert44     = require('gl-mat4/invert')
var cross        = require('gl-vec3/cross')
var normalize3   = require('gl-vec3/normalize')
var dot3         = require('gl-vec3/dot')

function len3(x, y, z) {
  return Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2))
}

function clamp1(x) {
  return Math.min(1.0, Math.max(-1.0, x))
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
    this.computedMatrix[i] = 0.5
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
  //Recompute curves
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
  var theta  = this.computedAngle[0]
  var phi    = this.computedAngle[1]

  var ctheta = Math.cos(theta)
  var stheta = Math.sin(theta)
  var cphi   = Math.cos(phi)
  var sphi   = Math.sin(phi)

  var center = this.computedCenter

  var wx = ctheta * cphi 
  var wy = stheta * cphi
  var wz = sphi

  var sx = -ctheta * sphi
  var sy = -stheta * sphi
  var sz = cphi

  var eye = this.computedEye
  var mat = this.computedMatrix
  for(var i=0; i<3; ++i) {
    var x      = wx * right[i] + wy * toward[i] + wz * up[i]
    eye[i]     = center[i] + radius * x
    mat[4*i+1] = sx * right[i] + sy * toward[i] + sz * up[i]
    mat[4*i+2] = x
    mat[4*i+3] = 0.0
  }

  var ax = mat[1]
  var ay = mat[5]
  var az = mat[9]
  var bx = mat[2]
  var by = mat[6]
  var bz = mat[10]
  var cx = ay * bz - az * by
  var cy = az * bx - ax * bz
  var cz = ax * by - ay * bx
  var cl = len3(cx, cy, cz)
  cx /= cl
  cy /= cl
  cz /= cl
  mat[0] = cx
  mat[4] = cy
  mat[8] = cz

  for(var i=0; i<3; ++i) {
    var rr = 0.0
    for(var j=0; j<3; ++j) {
      rr += mat[i+4*j] * eye[j]
    }
    mat[12+i] = -rr
  }
  mat[15] = 1.0
}

proto.tick = function(t) {
  this._recalcMatrix(t)
  this._lastTick = t
  this._isDirty = true
}

var SCRATCH0 = [0.1,0.1,0.1]
var SCRATCH1 = SCRATCH0.slice()
var SCRATCH2 = SCRATCH0.slice()
var SCRATCH_M = [0.1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

proto.setMatrix = function(t, matrix) {
  if(t < this._lastTick) {
    return
  }

  invert44(SCRATCH_M, matrix)
  var ex = SCRATCH_M[12]
  var ey = SCRATCH_M[13]
  var ez = SCRATCH_M[14]

  var dx = mat[2]
  var dy = mat[6]
  var dz = mat[10]
  var dl = len3(dx, dy, dz)
  dx /= dl
  dy /= dl
  dz /= dl

  this.radius.curve(t)
  var r  = Math.exp(this.computedRadius[0])
  var cx = ex + r * dx
  var cy = ey + r * dy
  var cz = ez + r * dz

  SCRATCH0[0] = cx
  SCRATCH0[1] = cy
  SCRATCH0[2] = cx

  SCRATCH1[0] = ex
  SCRATCH1[1] = ey
  SCRATCH1[2] = ez

  SCRATCH2[0] = mat[1]
  SCRATCH2[1] = mat[5]
  SCRATCH2[2] = mat[9]

  this.lookAt(t, SCRATCH0, SCRATCH1, SCRATCH2)
}

proto.rotate = function(t, dtheta, dphi) {
  this.radius.curve(t)
  var radius = Math.exp(this.computedRadius[0])
  this.radius.curve(this._lastTick)
  this.angle.move(t, dtheta / radius, dphi / radius)
}

proto.zoom = function(t, dz) {
  this.radius.move(t, dz)
}

proto.pan = function(t, dx, dy) {
  this._recalcMatrix(t)
  var mat = this.computedMatrix
  var rad = Math.exp(this.computedRadius[0])

  var ux = mat[1]
  var uy = mat[5]
  var uz = mat[9]
  var ul = len3(ux, uy, uz)
  ux /= ul
  uy /= ul
  uz /= ul

  var rx = mat[0]
  var ry = mat[4]
  var rz = mat[8]
  var ru = rx * ux + ry * uy + rz * uz
  rx -= ux * ru
  ry -= uy * ru
  rz -= uz * ru
  var rl = len3(rx, ry, rz)
  rx /= rl
  ry /= rl
  rz /= rl

  var vx = rx * dx + ux * dy
  var vy = ry * dx + uy * dy
  var vz = rz * dx + uz * dy
  this.center.move(t, rad * vx, rad * vy, rad * vz)

  this._recalcMatrix(this._lastTick)
}

//Recenters the coordinate axes
proto.tare = function(t, axes) {
  if(t < this._lastTick) {
    return
  }

  //Get the axes for tare
  var ushift, vshift
  if(typeof axes === 'number') {
    ushift = (axes)|0
    vshift = (ushift + 1) % 3
  } else if(typeof axes === 'string') {
    var laxes = axes.toUpperCase()
    ushift = laxes.charCodeAt(0) - 88
    if(axes.length === 2) {
      vshift = laxes.charCodeAt(1) - 88
    } else {
      vshift = (ushift + 1) % 3
    }
  } else {
    ushift = 1
    vshift = 0
  }
  if(ushift === vshift || 
    ushift < 0 || ushift >= 3 ||
    vshift < 0 || vshift >= 3) {
    return
  }

  //Recompute state for new t value
  this._recalcMatrix(t)

  //Read out matrix components
  var mat = this.computedMatrix

  //Get right and up vectors
  var ux = mat[ushift]
  var uy = mat[ushift+4]
  var uz = mat[ushift+8]
  var ul = len3(ux, uy, uz)
  ux /= ul
  uy /= ul
  uz /= ul

  var rx = mat[vshift]
  var ry = mat[vshift+4]
  var rz = mat[vshift+8]
  var ru = rx * ux + ry * uy + rz * uz
  rx -= ux * ru
  ry -= uy * ru
  rz -= uz * ru
  var rl = len3(rx, ry, rz)
  rx /= rl
  ry /= rl
  rz /= rl

  var eye    = this.computedEye
  var center = this.computedCenter
  var tx = eye[0] - center[0]
  var ty = eye[1] - center[1]
  var tz = eye[2] - center[2]
  var tl = len3(tx, ty, tz)
  var eu = (tx * ux + ty * uy + tz * uz) / tl
  var er = (tx * rx + ty * ry + tz * rz) / tl
  var ef = (tx * (uy * rz - uz * ry) +
            ty * (uz * rx - ux * rz) +
            tz * (ux * ry - uy * rx)) / tl

  this.radius.idle(t)
  this.center.idle(t)
  this.up.jump(t, ux, uy, uz)
  this.right.jump(t, rx, ry, rz)

  var phi, theta
  if(ushift === 2) {
    var cx = mat[1]
    var cy = mat[5]
    var cz = mat[9]
    var cr = (cx * rx + cy * ry + cz * rz)
    var cf =-(cx * (uy * rz - uz * ry) +
              cy * (uz * rx - ux * rz) +
              cz * (ux * ry - uy * rx))
    if(eu < 0) {
      phi = -Math.PI/2
    } else {
      phi = Math.PI/2
    }
    theta = Math.atan2(cf, cr)
  } else {
    phi = Math.asin(clamp1(eu))
    theta = Math.atan2(ef, er)
  }

  this.angle.jump(t, theta, phi)

  //Reset state of coordinates
  this._recalcMatrix(this._lastTick)
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

proto.lookAt = function(t, eye, center, up) {
  if(t < this._lastTick) {
    return
  }
  this._recalcMatrix(t)

  eye    = eye    || this.computedEye
  center = center || this.computedCenter
  up     = up     || this.computedUp
  
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
  this.radius.set(t, Math.log(tl))

  var fx = uy * rz - uz * ry
  var fy = uz * rx - ux * rz
  var fz = ux * ry - uy * rx
  var fl = len3(fx, fy, fz)
  fx /= fl
  fy /= fl
  fz /= fl

  var tu = ux*tx + uy*ty + uz*tz
  var tr = rx*tx + ry*ty + rz*tz
  var tf = fx*tx + fy*ty + fz*tz

  var phi   = Math.asin(clamp1(tu))
  var theta = Math.atan2(tf, tr)

  var angleState = this.angle._state
  var lastTheta  = angleState[angleState.length-1]
  var lastPhi    = angleState[angleState.length-2]
  lastTheta      = lastTheta % (2.0 * Math.PI)
  var dp = Math.abs(lastTheta + 2.0 * Math.PI - theta)
  var d0 = Math.abs(lastTheta - theta)
  var dn = Math.abs(lastTheta - 2.0 * Math.PI - theta)
  if(dp < d0) {
    lastTheta += 2.0 * Math.PI
  }
  if(dn < d0) {
    lastTheta -= 2.0 * Math.PI
  }

  this.angle.jump(this.angle.lastT(), lastTheta, lastPhi)
  this.angle.set(t, theta, phi)

  this._recalcMatrix(this._lastTick)
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
    phi    = Math.acos(ut)
    theta  = Math.acos(rt)
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