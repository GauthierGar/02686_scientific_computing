function [xdot,J] = CSTR3d(t,x,u,p)

% x = T
% u = F
% p = [CAin; CBin; Tin;  ....] i.e all the parameters

T = x(3);
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

CA = x(1);%CAin+(1/beta)*(Tin-T); % Equation (21a) in ACODSpaper
CB = x(2);%CBin+(2/beta)*(Tin-T); % Equation (21b) in ACODSpaper

k = k0*exp(-EaR*(1/T)); % Equation (3) ACODSpaper
r = k*CA*CB;  % Equation (2) ACODSpaper
RT = beta*r; % Equation (5) ACODSpaper
RA = -r; % Equation 5a ECC2020
RB = -2*r; %Equation 5b ECC2020

Tdot = zeros(3,1);
Tdot(1) = (F/V)*(CAin - CA) + RA; % Equation 6a ECC2020
Tdot(2) = (F/V)*(CBin - CB) + RB; % Equation 6b ECC2020
Tdot(3) = (F/V)*(Tin - T) + RT; %

xdot = Tdot;

J= zeros(3,3);
J(1,1)=-(F/V)-k*CB;
J(2,1)=-k*CA;
J(3,1)=-EaR*r*(1/(T.^2));
J(1,2)=-2*k*CB;
J(2,2)=-(F/V)-2*k*CA;
J(3,2)=-2*EaR*r*(1/(T.^2));
J(1,3)=beta*k*CB;
J(2,3)=beta*k*CA;
J(3,3)=-(F/V)+EaR*beta*r*(1/(T.^2));

end