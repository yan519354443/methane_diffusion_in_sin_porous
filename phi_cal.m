function result = phi_cal(rho, T, Tc,Gc)
    result = zeros(size(rho));
    [p,m,n] = size(rho);
    Tr = T/Tc;
    for i = 1:m
        for j = 1:n
            result(1,i,j) = 2.0*rho(1,i,j)*((8*Tr/(3-rho(1,i,j))-3*rho(1,i,j))/8/Tr-1/3.0)/6.0/Gc;
        end
    end
end