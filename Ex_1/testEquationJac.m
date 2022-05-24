function [f,J] = testEquationJac(t,x,varargin)
% The test function
%f = exp(-t);
f = -x;
% The Jaccobian
J = x;
end