function [Q] = midpoint(P)
% MIDPOINT finds a vector, Q, with elements equal to the mid pointes of successsive values of P. 
Q = (P(2:end)-P(1:end-1))./2;
