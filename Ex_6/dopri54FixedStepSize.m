function [T,X] = dopri54FixedStepSize(func,t0,tN,N,x0,varargin)
    % func  :   The function to evaluate
    % t0,tN :   The initial and final time
    % N     :   Number of steps
    % x0    :   The initial value

    % Define the variables
    h = (tN-t0)/N;
    nx = size(x0,1);
    X = zeros(nx,N+1);
    T = zeros(1,N+1);
    
    % Begin DOPRI5(4) algorithm
    X(:,1) = x0;
    T(:,1) = t0;

    for k=1:N
        [tf,xf] = dopri54Step(func,T(:,k),X(:,k),h,varargin{:});
        % We add them in the list
        T(:,k+1) = tf;
        X(:,k+1) = xf;
    end 
    
    % We return a good form for results 
    T=T';
    X=X';
end 