function xdot = CSTR1d(t,x,u,p)

% x = Temperature
% u = Flow
% p = [CAin; CBin; Tin;  ....] i.e all the parameters

T = x; % The temperature
F = u; % The flow rate 
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

CA = CAin+(1/beta)*(Tin-T); % Equation (21a) in ACODS2020
CB = CBin+(2/beta)*(Tin-T); % Equation (21b) in ACODS2020

k = k0*exp(-EaR*(1/T)); % Equation (3) ACODS2020
r = k*CA*CB;  % Equation (2) ACODS2020
RT = beta*r; % Equation (5) ACODS2020

tdot = (F/V)*(Tin - T) + RT; % Equation (20b) in ACODS2020
xdot = tdot;
end