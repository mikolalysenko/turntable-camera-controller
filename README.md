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
* `options.up` the up vector of the camera
* `options.right` the right vector the camera
* `options.radius` the distance from the camera to the objective
* `options.theta` the longitudinal angle of the camera
* `options.phi` the latitudinal angle of the camera

**Returns** A new camera controller

## Properties

#### `controller.center`

#### `controller.up`

#### `controller.right`

#### `controller.radius`

#### `controller.angle`


## 3d-camera-core interface

#### `controller.dirty()`

#### `controller.get(matrix)`


## Interaction

#### `controller.translate(t, x, y, z)`

#### `controller.pan(t, x, y)`

#### `controller.zoom(t, r)`

#### `controller.rotate(t, dx, dy)`

# License
(c) 2015 Mikola Lysenko. MIT License