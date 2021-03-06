function [T,X,info] = dopri54AdaptativeStepSize(func,t0,tf,x0,N,abstol,reltol,varargin)
    
    %Error controller parameters
    epstol = 0.8;   % target
    facmin = 0.1;   % maximum decrease factor
    facmax = 5.0;   % maximum increase factor
    
    %Initial conditions
    t = t0;
    h = (tf-t0)/N;
    x = x0;

    %Output
    T=t;
    X=x';
    
    % Information output
    nfun = 0;
    
    %% Main algo
    while t < tf
        % The last step
        if (t+h>tf)
            h = tf-t;
        end

        AcceptStep = false;
        while ~AcceptStep
            %Take step of size h
            [~,x1] = dopri54Step(func,t,x,h,varargin{:});
            nfun = nfun + 7;

            %Take step of size h/2
            hm = 0.5*h;
            [tm,xm] = dopri54Step(func,t,x,hm,varargin{:});
            nfun = nfun + 7;
            [~,x1hat] = dopri54Step(func,tm,xm,hm,varargin{:});
            nfun = nfun + 7;

            %Error estimation
            e = x1hat-x1;
            r = max(abs(e)./max(abstol,abs(x1hat).*reltol));
            
            % If we accept the step then
            AcceptStep = (r <= 1.0);
            if AcceptStep
                t = t+h;
                x = x1hat;

                T = [T;t];
                X = [X;x'];        
            end
            % Asymptotic step size controller
            h = max(facmin,min(sqrt(epstol/r),facmax))*h;
        end
    end
    info.nfun = nfun;
end 