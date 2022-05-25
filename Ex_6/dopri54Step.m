function [tf, xf] = dopri54Step(func,t,x,h,varargin)
    
    % The variable of the DOPRI5(4)
    a2 = 1/5;
    a3 = [3/40 9/40];
    a4 = [44/45 -56/15 32/9];
    a5 = [19872/6561 -25360/2187 64448/6561 -212/729];
    a6 = [9017/3168 -355/33 46732/5247 49/176 -5103/18656];
    a7 = [35/384 0 500/1113 125/192 -2187/6784 11/84];
    
    % The different operations of DOPRI algorithm
    x1 = feval(func,t,x,varargin{:});
    x2 = feval(func,t+h/5,x+h*(x1*a2),varargin{:});
    x3 = feval(func,t+(3*h/10),x+h*(x1*a3(1)+x2*a3(2)),varargin{:});
    x4 = feval(func,t+(4*h/5),x+h*(x1*a4(1)+x2*a4(2)+x3*a4(3)),varargin{:});
    x5 = feval(func,t+(8*h/9),x+h*(x1*a5(1)+x2*a5(2)+x3*a5(3)+x4*a5(4)),varargin{:});
    x6 = feval(func,t+h,x+h*(x1*a6(1)+x2*a6(2)+x3*a6(3)+x4*a6(4)+x5*a6(5)),varargin{:});
    
    % The final value of x and t
    tf = t+h;
    xf = x+h*(a7(1)*x1+a7(2)*x2+a7(3)*x3+a7(4)*x4+a7(5)*x5+a7(6)*x6);
end