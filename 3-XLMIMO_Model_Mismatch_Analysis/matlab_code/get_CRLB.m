function c = get_CRLB(c)
    if c.wave_type == "PWM"
        c = get_CRLB_PWM(c);
    elseif c.wave_type == "SWM"
        c = get_CRLB_SWM(c);
    elseif c.wave_type == "SWM_A1"
        c = get_CRLB_SWM_A1(c);
    end

end