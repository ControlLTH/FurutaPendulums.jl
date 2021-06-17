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
    d_current::AnalogIn
    q_current::AnalogIn
    torque::AnalogOut # Uses internal FOC controller which does not work well
    voltage::AnalogOut

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
            AnalogIn(m, Unsigned(47)),
            AnalogIn(m, Unsigned(48)),
            AnalogOut(m, Unsigned(40)),
            AnalogOut(m, Unsigned(41)),
        )
    end
end
