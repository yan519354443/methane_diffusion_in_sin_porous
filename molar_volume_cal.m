function result = molar_volume_cal(pressure, R, a, b, T)  
% https://www.youtube.com/watch?v=2pnTsHo2qSM
    result = zeros(size(pressure));
    [p,m,n] = size(pressure);
    initial_guess = 3e-3;
    for i = 1:m
        for j = 1:n
            pp = pressure(1,i,j);   %pp in mpa
            func = @(v)R*T*v.*(v+b).*sqrt(T)-a*(v-b)-pp*(v-b).*v.*(v+b).*sqrt(T);
            result(1,i,j) = fzero(func, initial_guess)*10e-9; 
        end
    end
end