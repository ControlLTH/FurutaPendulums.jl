# FurutaPendulums

This package implements the [AbstractControlProcesses](https://github.com/ControlLTH/AbstractControlProcesses.jl) interface to provide connections to both the physical furuta pendulum as well as a simulated version of it.

## Installation
First you need moberg installed and configure on your computer.

To install the package after that you need to run
```julia
import Pkg
Pkg.add("https://github.com/ControlLTH/FurutaPendulums.jl")
```

## Usage

```julia
using FurutaPendulums
furuta = FurutaPendulum() # or SimulatedFurutaPendulum()
initialize(furuta)
# Initialize time for periodic stepping to work both with simulated and physical process
last_time = periodic_wait(furuta, 0, 0) 

K = ...

dt = 0.01
for i in 1:1000
    x = measure(furuta)
    # Calculate control signal
    u = -K*x
    control(furuta, u)

    last_time = periodic_wait(furuta, last_time, dt)
end
```
