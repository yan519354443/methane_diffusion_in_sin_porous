function result = knusen_diffusion(rho, R, T, molecular_weight, pore_width)
    result = zeros(size(rho));
    [p,m,n] = size(rho);
    for i = 1:m
        for j = 1:n
            result(1,i,j) = (pore_width(i,1)/3)*sqrt(8*R*T/(pi*molecular_weight));
        end
    end
end