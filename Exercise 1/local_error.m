function [l] = local_error(mesh,approx)
    n = length(mesh);
    
    exa = zeros(n,1);
    for i=1:(n-1)
        exa(i) = exact_local(approx(i+1),mesh(i+1),mesh(i));
    end
    
    l = approx(2:n) - exa;
end