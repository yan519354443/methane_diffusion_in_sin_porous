function [result1, result2, result3] = eos_parameters_cal(Tc, Pc, Vm_cr, Zcr)

    result1 = Pc * Vm_cr / (Zcr * Tc);  % R
    result2 = 0.42748 * (result1^2) * (Tc^2.5) / Pc; % a
    result3 = 0.08644 * result1 * Tc / Pc; % b

end