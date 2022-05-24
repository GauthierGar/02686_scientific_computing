function [T,X] = ExplicitEulerFixedStepSize(fun,t0,tN,N,x0,varargin)

% Define the variables
dt = (tN-t0)/N;
nx = size(x0,1);
X = zeros(nx,N+1);
T = zeros(1,N+1);

% Euler Explicit Method 
T(:,1) = t0;
X(:,1) = x0;

for k=1:N
    f = feval(fun,T(:,k),X(:,k),varargin{:});
    T(:,k+1) = T(:,k) + dt;
    X(:,k+1) = X(:,k) + f*dt;
end 

% Return a good form for result 
T = T';
X = X';
end

