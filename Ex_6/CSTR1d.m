function [xdot,J] = CSTR1d(t,x,u,p)

% x = T
% u = F
% p = [CAin; CBin; Tin;  ....] i.e all the parameters

T = x;
F = u;
rho = p(1);
cp = p(2);
k0 = p(3);
EaR = p(4);
delHr = p(5);
V = p(6);
CAin = p(7); % you receive that from the parameter vector, p
CBin = p(8);  % you receive that from the parameter vector, p
Tin = p(9); % you receive that from the parameter vector, p

beta = -delHr/(rho*cp);

CA = CAin+(1/beta)*(Tin-T); % Equation (21a) in ACODSpaper
CB = CBin+(2/beta)*(Tin-T); % Equation (21b) in ACODSpaper

k = k0*exp(-EaR*(1/T)); % Equation (3) ACODSpaper
r = k*CA*CB;  % Equation (2) ACODSpaper
RT = beta*r; % Equation (5) ACODSpaper

Tdot = (F/V)*(Tin - T) + RT; %
xdot = Tdot;

J= zeros(1,1);
J(1,1)=-(F/V)-k*CB;

end