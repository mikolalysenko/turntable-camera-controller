# turntable-camera-controller

A low level controller for a turntable camera with input interpolation.  This module is compatible with the controller interface for 3d-camera-core.

# Example

```javascript
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

#### `controller.get(matrix)`


## Interaction

#### `controller.zoom(t, dr)`

#### `controller.pan(t, dx, dy)`

#### `controller.rotate(t, dx, dy)`



#### `controller.setEye(t, x, y, z)`

#### `controller.moveEye(t, dx, dy, dz)`

#### `controller.lookAt(t, center, up, eye)`

#### `controller.setUp(t, up)`

#### `controller.recenter(t)`


## State logging

#### `controller.flush(t)`

#### `controller.idle(t)`

# License
(c) 2015 Mikola Lysenko. MIT License