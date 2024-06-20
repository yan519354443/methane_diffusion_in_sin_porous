function result = pressure_cal(rho, phi, Gc)
    result = zeros(size(rho));
    [p,m,n] = size(rho);
    for i = 1:m
        for j = 1:n
        result(1,i,j) = rho(1,i,j)/3 + 3*Gc*(phi(1,i,j)^2);
        end
    end
end