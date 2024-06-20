function result = transition_diffusion(md, kd)
    result = zeros(size(md));
    [p,m,n] = size(md);
    for i = 1:m
        for j = 1:n
            result(1,j,i) = (1/md(1,i,j) + 1/kd(1,i,j))^-1;
        end
    end
end