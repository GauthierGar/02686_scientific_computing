function [x,nfunNewt] = NewtonsMethodODEAdaptative(func, tk, xk, dt, xinit, tol, maxit, varargin)
    
    k = 0;
    t = tk + dt;
    x = xinit;
    [f,J] = feval(func,t,x,varargin{:});
    R = x - f*dt - xk;
    I = eye(length(xk));

    nfunNewt = 0;
    
    while( (k<maxit) && (norm(R,'inf')>tol) )
        k = k+1;
        dRdx = I - J*dt;
        dx = dRdx\R;
        x = x - dx;
        [f,J] = feval(func,t,x,varargin{:});
        R = x - dt*f - xk;
        
        nfunNewt = nfunNewt + 1;
    end
    
end 
