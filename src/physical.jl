import MobergIO: Moberg, DigitalOut, DigitalIn, AnalogOut, AnalogIn, EncoderIn
 
struct FurutaPendulum <: PhysicalProcess
    # Connections
    energy_limited::DigitalIn
    calibrate::DigitalOut
    system_reset::DigitalOut
    base_encoder::EncoderIn
    arm_encoder::EncoderIn
    base_angle::AnalogIn
    arm_angle::AnalogIn
    base_velocity::AnalogIn
    arm_velocity::AnalogIn
    total_energy::AnalogIn
    base_energy::AnalogIn
    arm_energy::AnalogIn
    control_signal::AnalogOut

    # Other data
    h::Float64

    function FurutaPendulum()
        m = Moberg()
        return new(
            DigitalIn(m, Unsigned(40)),
            DigitalOut(m, Unsigned(40)),
            DigitalOut(m, Unsigned(41)),
            EncoderIn(m, Unsigned(40)),
            EncoderIn(m, Unsigned(41)),
            AnalogIn(m, Unsigned(40)),
            AnalogIn(m, Unsigned(41)),
            AnalogIn(m, Unsigned(42)),
            AnalogIn(m, Unsigned(43)),
            AnalogIn(m, Unsigned(44)),
            AnalogIn(m, Unsigned(45)),
            AnalogIn(m, Unsigned(46)),
            AnalogOut(m, Unsigned(40)),
            0.001 # TODO look up correct sampletime
        )
    end
end


Base.show(io::IO, furuta::FurutaPendulum) = print(io, "FurutaPendulum()") # Should maybe check connection and print <Connected>/<Disconnected>

# calibrate!(p::FurutaPendulum) = write(p.calibrate, true)
# reset!(p::FurutaPendulum) = write(p.system_reset, true)

# set_torque!(p::FurutaPendulum, torque) = write(p.torque, torque)

# get_energy_flag(p::FurutaPendulum) = read(p.energy_limited)
# get_total_energy(p::FurutaPendulum) = read(p.total_energy)

# get_phi(p::FurutaPendulum) = read(p.base_angle)
# get_dphi(p::FurutaPendulum) = read(p.base_velocity)
# get_theta(p::FurutaPendulum) = read(p.arm_angle)
# get_dtheta(p::FurutaPendulum) = read(p.arm_velocity)

# get_phi_count(p::FurutaPendulum) = read(p.base_encoder) # 20000 counts / 2pi
# get_theta_count(p::FurutaPendulum) = read(p.arm_encoder) # 16384 counts / 2pi 
# 
# get_x(p::FurutaPendulum) = (get_phi(p), get_dphi(p), get_theta(p), get_dtheta(p))