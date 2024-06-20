function result = molecular_diffution(number_density, kb,mass_molecule, T, d)
    result = zeros(size(number_density));
    [p,m,n] = size(number_density);
    for i = 1:m
        for j = 1:n
            n_star = number_density(1,i,j)*(d^3);
            result(1,i,j) = 3/(8*number_density(1,i,j)*(d^2))*sqrt(kb*T/(pi*mass_molecule))*...
                (1 - n_star/1.09)*(1 + (n_star^2)*(0.4 - 0.83*(n_star^2)));
        end
    end
end