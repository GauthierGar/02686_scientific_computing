function xdot = CSTR3d(t,x,u,p)

% x = T
% u = F
% p = [CAin; CBin; Tin;  ....] i.e all the parameters

F = u; % The flow
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

CA = x(1); % Concentration of CA
CB = x(2); % Concentration of CB
T = x(3); % The temperature 

k = k0*exp(-EaR*(1/T));
r = k*CA*CB;
RT = beta*r;
RA = -r; 
RB = -2*r;

Tdot = zeros(3,1);
Tdot(1) = (F/V)*(CAin - CA) + RA; % Equation 6a ECC2020
Tdot(2) = (F/V)*(CBin - CB) + RB; % Equation 6b ECC2020
Tdot(3) = (F/V)*(Tin - T) + RT; % Equation (20b) in ACODS2020

xdot = Tdot;
end