'use strict'

var now            = require('right-now')
var bunny          = require('bunny')
var perspective    = require('gl-mat4/perspective')
var fit            = require('canvas-fit')
var createContext  = require('gl-context')
var createAxes     = require('gl-axes')
var createMesh     = require('gl-simplicial-complex')
var createCamera   = require('../turntable')

//Set up WebGL
var canvas = document.createElement('canvas')
document.body.appendChild(canvas)
window.addEventListener('resize', fit(canvas), false)
var gl = createContext(canvas, {}, render)

//Create objects for rendering
var bounds = [[-10,-10,-10], [10,10,10]]
var mesh = createMesh(gl, {
    cells: bunny.cells,
    positions: bunny.positions,
    colormap: 'jet'
  })
var axes = createAxes(gl, {
    bounds: bounds,
    tickSpacing: [1,1,1],
    textSize: 0.05
  })

//Set up camera
var projectionMatrix = new Array(16)
var camera = createCamera({
  radius: 50,
  center:  [
    0.5*(bounds[0][0]+bounds[1][0]),
    0.5*(bounds[0][1]+bounds[1][1]),
    0.5*(bounds[0][2]+bounds[1][2]) ]
})

var lastX = 0, lastY = 0
canvas.addEventListener('mousemove', function(ev) {
  if(ev.which) {
    camera.rotate(now(), 
      -150 * (ev.x - lastX) / gl.drawingBufferWidth, 
      150 * (ev.y - lastY) / gl.drawingBufferHeight)
  }
  lastX = ev.x
  lastY = ev.y
  console.log(ev)
})

//Redraw frame
function render() {
  //Update camera parameters
  var t = now()
  var delay = 30
  camera.idle(t - delay)
  camera.tick(t - 2 * delay)
  camera.flush(t - 10 * delay)

  //Compute parameters
  var cameraParams = {
    view: camera.get(),
    projection: perspective(
      projectionMatrix, 
      Math.PI/4.0,
      gl.drawingBufferWidth/gl.drawingBufferHeight,
      0.1,
      1000.0)
  }

  //Draw everything
  gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight)
  gl.enable(gl.DEPTH_TEST)
  axes.draw(cameraParams)
  mesh.draw(cameraParams)
}