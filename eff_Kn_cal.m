function result = eff_Kn_cal(rho, mass_molecule, d, pore_width)
    [p,m,n] = size(rho);
    Kn = zeros(size(rho));
    result = zeros(size(rho));
    for i = 1:m
        for j = 1:n
            Kn(1,i,j) = mass_molecule / (sqrt(2) * pi * rho(1,i,j) * (d^2)*pore_width(i));
            result(1,i,j) = Kn(1,i,j)*2 / pi * atan(sqrt(2) * Kn(1,i, j) ^ (-3 / 4));
        end
    end
end 