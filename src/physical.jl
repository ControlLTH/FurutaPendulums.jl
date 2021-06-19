import MobergIO: Moberg, DigitalOut, DigitalIn, AnalogOut, AnalogIn, EncoderIn
  
Base.@kwdef struct FurutaPendulum <: PhysicalProcess
    m::Moberg                   = Moberg()
    energy_limited::DigitalIn   = DigitalIn(m, Unsigned(40))
    calibrate::DigitalOut       = DigitalOut(m, Unsigned(40))
    system_reset::DigitalOut    = DigitalOut(m, Unsigned(41))
    base_encoder::EncoderIn     = EncoderIn(m, Unsigned(40))
    arm_encoder::EncoderIn      = EncoderIn(m, Unsigned(41))
    base_angle::AnalogIn        = AnalogIn(m, Unsigned(40))
    arm_angle::AnalogIn         = AnalogIn(m, Unsigned(41))
    base_velocity::AnalogIn     = AnalogIn(m, Unsigned(42))
    arm_velocity::AnalogIn      = AnalogIn(m, Unsigned(43))
    total_energy::AnalogIn      = AnalogIn(m, Unsigned(44))
    base_energy::AnalogIn       = AnalogIn(m, Unsigned(45))
    arm_energy::AnalogIn        = AnalogIn(m, Unsigned(46))
    d_current::AnalogIn         = AnalogIn(m, Unsigned(47))
    q_current::AnalogIn         = AnalogIn(m, Unsigned(48))
    torque::AnalogOut           = AnalogOut(m, Unsigned(40))  # Uses internal FOC controller which does not work well
    voltage::AnalogOut          = AnalogOut(m, Unsigned(41))

end
