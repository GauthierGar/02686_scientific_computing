function [T,X] = ImplicitEulerFixedStepSize(funJac,t0,tN,N,x0,varargin)

% Define the variables
dt = (tN-t0)/N;
nx = size(x0,1);
X = zeros(nx,N+1);
T = zeros(1,N+1);

tol = 1.0e-8;
maxit = 100;

% Euler Implicit Method 
T(:,1) = t0;
X(:,1) = x0;

for k=1:N
    f = feval(funJac,T(k),X(k),varargin{:});
    T(:,k+1) = T(:,k) + dt;
    xinit = X(:,k) + f*dt;
    X(:,k+1) = NewtonsMethodODE(funJac,T(:,k),X(:,k),dt,xinit,tol,maxit,varargin{:});
end 

% Return a good form for result 
T = T';
X = X';
end