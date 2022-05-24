function [e] = global_error(mesh,approx)
    e = exact(mesh) - approx;
end