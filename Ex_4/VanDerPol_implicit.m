function [xdot,Jac] = VanDerPol_implicit(t,x,mu)
     % VANDERPOL  Implementation of the Van der Pol model
     %
     % Syntax: xdot = VanDerPol(t,x,mu)
     xdot = zeros(2,1);
     xdot(1) = x(2);
     xdot(2) = (mu*(1-x(1)*x(1))*x(2))-x(1);
     
     Jac = zeros(2,2);
     Jac(2,1) = -2*mu*x(1)*x(2)-1.0;
     Jac(1,2) = 1.0;
     Jac(2,2) = mu*(1-x(1)*x(1));
end
