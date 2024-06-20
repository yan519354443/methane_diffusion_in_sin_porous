
T = 298; % in K
Tc = 190.4; % in K
Gc = -0.4; % the interaction strength of fluid particles
Pc = 4.595;   %Mpa
Vm_cr = 98.66; %m3/kg
Zcr = 0.287;
NA = 6.02*10^23;
kb = 1.3e-23; %Boltzmann constant
molecular_weight = 16.04; % g/mol
pressure = [[[0.331]]]; %MPa

[R, a , b] = eos_parameters_cal(Tc, Pc, Vm_cr, Zcr);  % R in J/mol*K
V = molar_volume_cal(pressure, R, a, b, T);  %v in m^3/kg