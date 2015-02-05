# turntable-camera-controller

A low level controller for a turntable camera with input interpolation.  This module is compatible with the controller interface for 3d-camera-core.

You probably don't need to use this directly in most cases.

# Example

```javascript
var now            = require('right-now')
var bunny          = require('bunny')
var perspective    = require('gl-mat4/perspective')
var fit            = require('canvas-fit')
var createContext  = require('gl-context')
var createAxes     = require('gl-axes')
var createMesh     = require('gl-simplicial-complex')
var createCamera   = require('turntable-camera-controller')

//Set up WebGL
var canvas = document.createElement('canvas')
document.body.appendChild(canvas)
window.addEventListener('resize', fit(canvas), false)
var gl = createContext(canvas, {}, render)

var controlDiv = document.createElement('div')
controlDiv.style.position = 'absolute'
controlDiv.style['z-index'] = 10
controlDiv.style.left = '10px'
controlDiv.style.top = '10px'
document.body.appendChild(controlDiv)

var delayControl = document.createElement('input')
delayControl.type = 'range'
delayControl.min = 0
delayControl.max = 200
delayControl.value = 30
controlDiv.appendChild(delayControl)

var tareButton = document.createElement('input')
tareButton.type = 'submit'
tareButton.value = 'Roll'
controlDiv.appendChild(tareButton)

var lookAtButton = document.createElement('input')
lookAtButton.type = 'submit'
lookAtButton.value = 'Reset'
controlDiv.appendChild(lookAtButton)


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

//Hook event listeners
var lastX = 0, lastY = 0

document.oncontextmenu = function(e) { 
  e.preventDefault()
  e.stopPropagation()
  return false 
}

canvas.addEventListener('mousemove', function(ev) {
  var dx = (ev.clientX - lastX) / gl.drawingBufferWidth
  var dy = (ev.clientY - lastY) / gl.drawingBufferHeight
  if(ev.which === 1) {
    camera.rotate(now(), 
      -150 * dx, 
       150 * dy)
  }
  if(ev.which === 3) {
    camera.pan(now(), -dx, dy)
  }
  lastX = ev.clientX
  lastY = ev.clientY
})

canvas.addEventListener('wheel', function(e) {
  camera.zoom(now(), e.deltaY)
})

tareButton.addEventListener('click', function() {
  camera.tare(now(), 2)
})

lookAtButton.addEventListener('click', function() {
  camera.lookAt(now(),
    [ 0, 0, -50 ],
    [ 0.5*(bounds[0][0]+bounds[1][0]),
      0.5*(bounds[0][1]+bounds[1][1]),
      0.5*(bounds[0][2]+bounds[1][2]) ],
    [ 0, 1, 0 ])
})

//Redraw frame
function render() {
  //Update camera parameters
  var t = now()
  var delay = +delayControl.value
  camera.idle(t - delay)
  camera.tick(t - 2 * delay)
  
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
```

# Install

```
npm i turntable-camera-controller
```

# API

## Constructor

#### `var controller = require('turntable-camera-controller')(options)`
Creates a new turntable controller with the given input parameters.

* `options.center` the center of the camera
* `options.eye` the location of the eye (optional)
* `options.up` the up vector of the camera (default [0,1,0])
* `options.right` the right vector the camera (default [1,0,0])
* `options.radius` the distance from the camera to the objective (default 1)
* `options.theta` the longitudinal angle of the camera (default
* `options.phi` the latitudinal angle of the camera

**Returns** A new camera controller

## 3d-camera-core interface

#### `controller.dirty()`
Test if the camera is dirty from last tick

#### `controller.get(matrix)`
Retrieves the current matrix for the camera controller.


## Interaction

#### `controller.idle(t)
Notifies the controller of being idle at time `t`

* `t` is the time at which input was idle

#### `controller.zoom(t, dr)`
Zooms in the camera at time `t`

* `t` is the time of the zoom event
* `dr` the change in zoom factor in log scale

#### `controller.pan(t, dx, dy)`
Pans the camera by the amount `dx` and `dy`

#### `controller.rotate(t, dx, dy)`
Rotates the camera

#### `controller.lookAt(t, center, up, eye)`
Similar semantics to `gl-mat4/lookAt`. 

* `t` is the time of the event
* `center` is the view target
* `up` is the upward axis (axis of rotation)
* `eye` is the position of the eye

#### `controller.tare(t[, axis])`
Moves the `up` vector of the camera to be relative to the current view

* `t` is the time of the event
* `axis` is an integer in [0,3) determining the axis of the new rotation system

#### `controller.setMatrix(t, matrix)`
Sets the camera controller state to match `matrix`

## State logging

#### `controller.flush(t)`
Clears out all events before time `t`

#### `controller.tick(t)`
Updates the camera to time `t`.

# License
(c) 2015 Mikola Lysenko. MIT License