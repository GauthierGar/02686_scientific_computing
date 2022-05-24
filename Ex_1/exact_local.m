function [a] = exact_local(approx,tk,tkmoins)
    a = exact(tk) - exact(tkmoins) + approx;
end
