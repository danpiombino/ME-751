function [q_dd,lam] = decz(z,nb)
% function [q_dd,lam] = decz(z)
% decomposes z into q_dd and lam
% 
% Inputs:   z: [r_dd;p_dd;lam]
%           nb: number of bodies
%
% Outputs:  q_dd: Accelerations in form [r_dd;p_dd]
%           lam: Lagrange multipliers
%
% Written By: Dan Piombino
% 11/3/16

q_dd = z(1:7*nb);
lam = z(7*nb+1:end);