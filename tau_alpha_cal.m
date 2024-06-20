function result = tau_alpha_cal(kne, md, kd, td, delta_t)
    result = zeros(size(md));
    [p,m,n] = size(md);
    for i = 1:m
        for j = 1:n
            if kne(1,i,j) <= 0.1
                result(1,i,j) = 3*md(1,i,j)/delta_t + 0.5;
            elseif kne(1,i,j) > 0.1 && kne(1,i,j) < 10
                result(1,i,j) = 3*td(1,i,j)/delta_t + 0.5;
            else
                result(1,i,j) = 3*kd(1,i,j)/delta_t + 0.5;
            end       
        end
    end

end