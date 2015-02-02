# turntable-controller

A low level controller for a turntable camera with input interpolation.

# Example

```javascript
```

# Install

```
npm i turntable-controller
```

# API

## Constructor

#### `var controller = require('turntable-controller')(options)`
Creates a new turntable controller with the given input parameters.

* `options.center` the center of the camera
* `options.up` the up vector of the camera
* `options.radius` the distance from the camera to the objective
* `options.theta` the longitudinal angle of the camera
* `options.phi` the latitudinal angle of the camera

**Returns** A new camera controller

## Properties

#### `controller.center`

#### `controller.up`

#### `controller.radius`

#### `controller.angle`


## Camera core interface

#### `controller.dirty()

#### `controller.get(matrix)`


## Interaction

#### `controller.translate(t, x, y, z)`

#### `controller.pan(t, x, y)`

#### `controller.zoom(t, r)`

#### `controller.rotate(t, dx, dy)`

# License
(c) 2015 Mikola Lysenko. MIT License