struct FurutaParams{T}
    α::T
    β::T
    γ::T
    δ::T
    τc::T
    τs::T
    h::T
    max_speed::T
    max_torque::T
    noise::T
end

mutable struct SimulatedFurutaPendulum{T, R<:AbstractRNG} <: SimulatedProcess
    params::FurutaParams{T}
    x::Vector{T} # θ, θdot, ϕ, ϕdot = base angle & speed, arm angle & speed
    u::T
    rng::R
end

"""
    SimulatedFurutaPendulum(; 
        T=Float64,
        J = T(1.54e-4),
        M = T(0),
        ma = T(0),
        mp = T(5.44e-3),
        la = T(4.3e-2),
        lp = T(6.46e-2),
        τc = T(0.0076),
        τs = T(0.008),
        h = T(0.002),
        max_speed = T(100),
        max_torque = T(1),
        rng = Random.GLOBAL_RNG,
    )

Creates an instance of a simulator for the furutapendulum. 
"""
function SimulatedFurutaPendulum(;
        T = Float64,
        J = T(1.54e-4),
        M = T(0),
        ma = T(0),
        mp = T(5.44e-3),
        la = T(4.3e-2),
        lp = T(6.46e-2),
        τc = T(0.0076),
        τs = T(0.008),
        h = T(0.001), # Simulation time, TODO find what this actually is
        max_speed = T(100),
        max_torque = T(1),
        x0 = zeros(T, 4),
        noise = T(0.01),
        rng = Random.GLOBAL_RNG,
    )
    g = 9.81
    α = J+(M+ma/3+mp)*la^2
    β = (M+mp/3)*lp^2
    γ = (M+mp/2)*la*lp
    δ = (M+mp/2)*g*lp
    
    params = FurutaParams{T}(α, β, γ, δ, τc, τs, h, max_speed, max_torque, noise)
    SimulatedFurutaPendulum{T, typeof(rng)}(params, x0, zero(T), rng)
end

function step!(p, dt)
    u = p.u
    x = p.x
    h = p.params.h 

    steps = round(Int, dt / h)
    steps ≈ dt / h || throw(ArgumentError("`dt` need to be a multiple of the simulation step `h`."))
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
end

# See https://portal.research.lu.se/portal/files/4453844/8727127.pdf 
# and https://ctms.engin.umich.edu/CTMS/index.php?example=MotorSpeed&section=SystemModeling
# for source of dynamic equations
function f(p::SimulatedFurutaPendulum, x, τ)
    # Base angle/speed, arm angle/speed, torque
    θ, θdot, ϕ, ϕdot = x

    @unpack α, β, γ, δ, τc, τs = p.params

    # Friction forces
    τF = if abs(x[4]) > 0.01 
        τc*sign(x[4]) 
    elseif abs(τ) < τs
        τ
    else
        τs*sign(τ)
    end

    cost = cos(θ)
    sint = sin(θ)

    ψ = 1 / (α*β - γ^2 + (β^2 + γ^2)*sint^2)
    dθ    = θdot
    dθdot = ψ * (β*(α+β*sint^2)*cost*sint*ϕdot^2
        + 2*β*γ*(1-sint^2)*sint*ϕdot*θdot - γ^2*cost*sint*θdot^2
        + δ*(α+β*sint^2)*sint - γ*cost*(τ-τF)) 
    dϕ    = ϕdot
    dϕdot = ψ*(β*γ*(sint^2-1)*sint*ϕdot^2 - 2*β^2*cost*sint*ϕdot*θdot
        + β*γ*sint*θdot^2 - γ*δ*cost*sint + β*(τ-τF))

    return [dθ, dθdot, dϕ, dϕdot]
end