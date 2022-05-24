function [f,J] = testEquationJac(t,x,varargin)
% The test function
f = exp(-t);
% The Jaccobian
J = -exp(-t);
end