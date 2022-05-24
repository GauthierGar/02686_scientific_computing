function [t1, x1] = ClassicalRungeKuttaStep(fun,t,x,h,varargin)

    h2 = 0.5*h;
    alpha = h/6;
    beta = h/3;
    
    T1 = t;
    X1 = x;
    [F1,J1] = feval(fun,T1,X1,varargin{:});

    T2 = x+h2;
    X2 = x+h2*F1;
    [F2,J2] = feval(fun,T2,X2,varargin{:});

    T3 = T2;
    X3 = x+h2*F2;
    [F3,J3] = feval(fun,T3,X3,varargin{:});

    T4 = t+h;
    X4 = x+h*F3;
    [F4,J4] = feval(fun,T4,X4,varargin{:});

    t1 = T4;
    x1 = x + alpha*(F1+F4) + beta*(F2+F3);
end