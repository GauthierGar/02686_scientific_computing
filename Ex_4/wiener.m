function [w] = wiener(n,n_w,dt)

w = zeros(n_w,n);

moyenne = 0;
sigma = dt;

w = normrnd(moyenne,sigma,n_w,n);

w = cumsum(w,2);

end