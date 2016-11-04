function [q_dd,lam] = decz(z,nb)
% function [q_dd,lam] = decz(z)
% decomposes z into q_dd and lam
% 
% Inputs:   z: [r_dd;p_dd;lam]
%           nb: number of bodies
%
% Outputs:  q_dd: Accelerations in form [r1_dd;p1_dd;r2_dd;p2_dd; ...]
%           lam: Lagrange multipliers
%
% Written By: Dan Piombino
% 11/3/16

nc = length(z)-8*nb;
r_dd = z(1:3*nb);
p_dd = z(3*nb+1:7*nb);
lam = z(7*nb+1:end);
q_dd = rp2q(r_dd,p_dd);