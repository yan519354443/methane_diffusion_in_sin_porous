function result = r_coefficient(kne, pore_width,width_real_each_grid)
    result = zeros(size(kne));
    [p,m,n] = size(kne);
    xigema = 1;
    a1 = (2 - xigema) / xigema * (1 - 0.1817 * xigema);
    a2 = 1 / pi + 1 / 2 * (a1 ^ 2); 

    for i = 1:m
        for j = 1:n
            result(1,i,j) = 1 / (1 + sqrt(pi / 6) * (((pore_width(i)/width_real_each_grid)) ^ 2 / (4 * kne(1,i,j))...
                                                + a1 + (2 * a2 - 8 / pi) * kne(1,i,j)));
        end
    end
end