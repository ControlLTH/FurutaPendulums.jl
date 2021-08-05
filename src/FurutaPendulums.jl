module FurutaPendulums

export FurutaPendulum, SimulatedFurutaPendulum, AbstractFurutaPendulum

using Random, UnPack

using Reexport                           # Use reexport to avoid writing all exports
@reexport using AbstractControlProcesses # Import and export all exported names from ACP
const ACP = AbstractControlProcesses     # Less to write, extending function with new

include("physical.jl")
include("simulated.jl")

const AbstractFurutaPendulum = Union{FurutaPendulum, SimulatedFurutaPendulum}

ACP.num_outputs(p::AbstractFurutaPendulum) = 4
ACP.num_inputs(p::AbstractFurutaPendulum) = 1
ACP.outputrange(p::AbstractFurutaPendulum) = [(0, 2pi), (-1e10, 1e10), (0, 2pi), (-1e10, 1e10)]
ACP.inputrange(p::AbstractFurutaPendulum) = [(-1, 1)]
ACP.isstable(p::AbstractFurutaPendulum) = false
ACP.isasstable(p::AbstractFurutaPendulum) = false
ACP.sampletime(p::FurutaPendulum) = 0.001 # TODO find the correct one
ACP.sampletime(p::SimulatedFurutaPendulum) = p.params.h
ACP.bias(p::AbstractFurutaPendulum) = 0

ACP.control(p::FurutaPendulum, u::Real) = write(p.voltage, u)
ACP.control(p::SimulatedFurutaPendulum, u::Real) = p.u = clamp(u, -p.params.max_torque, p.params.max_torque)
ACP.measure(p::FurutaPendulum) = [read(p.base_angle), read(p.base_velocity), read(p.arm_angle), read(p.arm_velocity)]
ACP.measure(p::SimulatedFurutaPendulum) = p.x .+ p.params.noise * randn(p.rng, 4)

function ACP.periodic_wait(p::SimulatedFurutaPendulum, last_time, dt)
    step!(p, dt)
    last_time + dt
end

# ACP.initialize(p::FurutaPendulum) = write(p.calibrate, true) 
ACP.initialize(p::SimulatedFurutaPendulum) = fill!(p.x, 0) # Used as a reset

end
