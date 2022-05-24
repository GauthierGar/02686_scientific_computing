function [T,X] = ClassicalRungeKuttaFixedStepSize(funJac,t0,tN,N,x0,varargin)

% Define the variables
dt = (tN-t0)/N;
nx = size(x0,1);
X = zeros(nx,N+1);
T = zeros(1,N+1);

% Euler Explicit Method 
T(:,1) = t0;
X(:,1) = x0;

for k=1:N
    [t,x] = ClassicalRungeKuttaStep(funJac,T(:,k),X(:,k),dt,varargin);
    T(:,k+1) = t;
    X(:,k+1) = x;
end

% Return a good form for result 
T = T';
X = X';
end