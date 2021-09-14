struct FurutaParams{T}
    α::T
    β::T
    γ::T
    δ::T
    τc::T
    τs::T
    τv::T
    h::T
    max_speed::T
    max_torque::T
    noise::T
    friction::Symbol
end

mutable struct SimulatedFurutaPendulum{T, R<:AbstractRNG} <: SimulatedProcess
    params::FurutaParams{T}
    x::Vector{T} # θ, θdot, ϕ, ϕdot = base angle & speed, arm angle & speed
    u::T
    rng::R
end

"""
    SimulatedFurutaPendulum(; 
        J = 1.54e-4,
        M = 0,
        ma = 0,
        mp = 5.44e-3,
        la = 4.3e-2,
        lp = 6.46e-2,
        friction = :none,
        τc = 0.0076,
        τs = 0.008,
        τv = 0.008,
        h = 0.002,
        x0 = zeros(4),
        max_speed = 100,
        max_torque = 1,
        noise = 0.01,
        friction = :viscous,
        rng = Random.GLOBAL_RNG,
    )

Creates an instance of a simulator for the furutapendulum. 
It is made up of two inertial bodies, a center pillar with 
intertia `J` connected to a horizontal arm with length `la`
and homogenuously distributed mass `ma`. The other end of the
arm is connected to the pendulum which has length `lp`, homogenuously 
distributed mass `mp` and point distributed mass `M` at the end.
The `ϕ`-joint is in the center pillar while the `θ`-joint is for the
pendulum arm. There is significant `friction` in the `ϕ`-joint which 
is modeled using either coulomb and viscous friction `:viscous`, 
or coulomb friction and stiction `:stiction` (for more info see Gäfvert
https://portal.research.lu.se/portal/files/4453844/8727127.pdf).
The simulation takes steps of `h` seconds, and the initial value `x0`
is for `[ϕ, ϕdot, θ, θdot]`.
There is `max_speed` and `max_torque` limiting the speed and torque, and
`noise` that adds white noise with this scaling factor on top of the 
measurements.
"""
function SimulatedFurutaPendulum(;
        J = 1.54e-4,
        M = 0,
        ma = 0,
        mp = 5.44e-3,
        la = 4.3e-2,
        lp = 6.46e-2,
        τc = 0.0076,
        τs = 0.008,
        τv = 0.008,
        h = 0.002, # Simulation time, TODO find what this actually is
        max_speed = 100,
        max_torque = 0.04,
        x0 = [0.0, 0.0, π, 0.0], # Start down as default
        noise = 0.01,
        friction::Symbol = :viscous,
        rng = Random.GLOBAL_RNG,
    )
    g = 9.81
    α = J+(M+ma/3+mp)*la^2
    β = (M+mp/3)*lp^2
    γ = (M+mp/2)*la*lp
    δ = (M+mp/2)*g*lp
    
    params = FurutaParams{Float64}(α, β, γ, δ, τc, τs, τv, h, max_speed, max_torque, noise, friction)
    SimulatedFurutaPendulum{Float64, typeof(rng)}(params, x0, 0.0, rng)
end

function step!(p::SimulatedFurutaPendulum, dt)
    u = p.u
    x = p.x
    h = p.params.h 

    steps = round(Int, dt / h)
    steps ≈ dt / h || throw(ArgumentError("`dt`=$(dt) need to be a multiple of the simulation step `h`=$(h)."))
    for _ in 1:steps
        #RK3/8
        k1 = f(p, x, u)
        k2 = f(p, x + h * k1 / 3, u)
        k3 = f(p, x + h * (-k1 / 3 + k2), u)
        k4 = f(p, x + h * (k1 - k2 + k3), u)
        x .+= h .* (k1 .+ 3 .* k2 .+ 3 .* k3 .+ k4) ./ 8
    end

    x[1] = rem2pi(x[1], RoundToZero)
    x[3] = rem2pi(x[3], RoundToZero)

    # Should this be here?
    x[2] = clamp(x[2] , -p.params.max_speed, p.params.max_speed)
    x[4] = clamp(x[4] , -p.params.max_speed, p.params.max_speed)

    nothing
end

# See https://portal.research.lu.se/portal/files/4453844/8727127.pdf 
# and https://ctms.engin.umich.edu/CTMS/index.php?example=MotorSpeed&section=SystemModeling
# for source of dynamic equations
function f(p::SimulatedFurutaPendulum, x, τ)
    # Base angle/speed, arm angle/speed, torque
    ϕ, ϕdot, θ, θdot = x

    @unpack α, β, γ, δ, friction = p.params

    # Friction forces
    τF = 0
    if friction === :viscous
        τF = viscous_friction(p, ϕdot)
    elseif friction === :stiction
        τF = stiction_friction(p, ϕdot, τ)
    end

    cost = cos(θ)
    sint = sin(θ)

    ψ = 1 / (α*β - γ^2 + (β^2 + γ^2)*sint^2)
    dϕ    = ϕdot
    dϕdot = ψ*(β*γ*(sint^2-1)*sint*ϕdot^2 - 2*β^2*cost*sint*ϕdot*θdot
        + β*γ*sint*θdot^2 - γ*δ*cost*sint + β*(τ-τF))
    dθ    = θdot
    dθdot = ψ * (β*(α+β*sint^2)*cost*sint*ϕdot^2
        + 2*β*γ*(1-sint^2)*sint*ϕdot*θdot - γ^2*cost*sint*θdot^2
        + δ*(α+β*sint^2)*sint - γ*cost*(τ-τF)) 

    return [dϕ, dϕdot, dθ, dθdot]
end

function viscous_friction(p::SimulatedFurutaPendulum, ϕdot)
    p.params.τc * sign(ϕdot) + p.params.τv * ϕdot
end
function stiction_friction(p::SimulatedFurutaPendulum, ϕdot, τ)
    if abs(ϕdot) > 0.02 
        p.params.τc*sign(ϕdot) 
    elseif abs(τ) < p.params.τs
        τ
    else
        τs*sign(τ)
    end
end

struct ValueContainer 
    v::Float64
end

function Base.getproperty(p::SimulatedFurutaPendulum, s::Symbol)
    K0 = 0.0000652 # kinetic energy constant base: 0.5*moment of inertia, kgm^2 
    K1 = 0.00000387 # kinetic energy constant pendulum 

    if s === :arm_energy
        ϕ, ϕdot, θ, θdot = p.x
        arm_energy = (K1 * θdot^2 + 0.0054 * 9.8 * 0.07/2 * (1 - cos(θ)))
        return ValueContainer(arm_energy)
    elseif s === :total_energy
        ϕ, ϕdot, θ, θdot = p.x
        arm_energy = (K1 * θdot^2 + 0.0054 * 9.8 * 0.07/2 * (1 - cos(θ)))
        base_energy = K0 * ϕdot^2
        return ValueContainer(base_energy + arm_energy)
    else
        return getfield(p, s)
    end
end

Base.read(c::ValueContainer) = c.v