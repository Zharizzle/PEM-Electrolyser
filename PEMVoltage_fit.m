function V_cell = PEMVoltage_fit(I, T, alpha_a, alpha_c, j0_a, j0_c, i_lim, R, F, t_mem, a)

    Vrev = 1.23;

    a_c = min(max(a, 0), 3);
    if a_c <= 1
        lambda = 0.043 + 17.81*a_c - 39.85*a_c^2 + 36.0*a_c^3;
    else
        lambda = 14 + 1.4*(a_c - 1);
    end
    sigma_raw = (0.005139*lambda - 0.00326) * exp(1268*(1/303 - 1/T));
    sigma = max(sigma_raw, 1e-6);

    eta_act  = 0;
    eta_ohm  = 0;
    eta_diff = 0;

    if I > 0
        eta_act = (R*T/(alpha_a*F)) * asinh(I/(2*j0_a)) + ...
                  (R*T/(alpha_c*F)) * asinh(I/(2*j0_c));
        eta_ohm = I * (t_mem/sigma);
        i_safe   = min(I, 0.999*i_lim);
        eta_diff = (R*T/F) * log(i_lim/(i_lim - i_safe));
    end

    V_cell = Vrev + eta_act + eta_ohm + eta_diff;
end
