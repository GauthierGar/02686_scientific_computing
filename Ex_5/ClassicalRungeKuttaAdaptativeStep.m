function [T,X,info] = ClassicalRungeKuttaAdaptativeStep(func,t0,tf,x0,N,abstol,reltol,varargin)

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
    hvec = [];
    rvec = [];
    err = [];
    rr = [];
    hh =[];
    
    %% Main algo
    while t < tf
        if (t+h>tf)
            h = tf-t;
        end

        AcceptStep = false;
        while ~AcceptStep
            %Take step of size h
            [~,x1] = ClassicalRungeKuttaStep(func,t,x,h,varargin{:});

            %Take step of size h/2
            hm = 0.5*h;
            [tm,xm] = ClassicalRungeKuttaStep(func,t,x,hm,varargin{:});
            [~,x1hat] = ClassicalRungeKuttaStep(func,tm,xm,hm,varargin{:});

            %Error estimation
            e = x1hat-x1;
            r = max(abs(e)./max(abstol,abs(x1hat).*reltol));

            AcceptStep = (r <= 1.0);
            if AcceptStep
                t = t+h;
                x = x1hat;

                T = [T;t];
                X = [X;x'];
                
                err = [err;e.'];
                rr = [rr;r.'];
                hh = [hh;h.'];

            end
            hvec(end+1) = h;
            rvec(end+1) = r;
            %asymptotic step size controller
            h = max(facmin,min(sqrt(epstol/r),facmax))*h;
        end
    end

    info.nfun = nfun;
    info.naccept = length(T);
    info.nreject = length(hvec) - length(T) + 1;
    info.nstep = length(hvec) + length(rvec);
    info.hvec = hvec;
    info.rvec = rvec;
    info.err = err;
    info.rr = rr;
    info.hh = hh;
end 