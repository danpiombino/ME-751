function [q_dd,lam] = maf(M,Jp,phi_q,F,tau,gamma)
% function [q_dd,lam] = maf(M,Jp,phi_q,F,tau,gamma)
% calculates "a = inv(m)*F"
% NOTE: Constraints MUST be in order of: all kinematic
%       constraints, then all euler parameterization constraints
% 
% Inputs:   M: Mass matrix [3*nb x 3*nb]
%           Jp: Euler Parameter Inertia Matrix [4*nb x 4*nb]
%           F: External Forces in LRF [3*nb x 1]
%           tau: External Torques in LRF [4*nb x 1]
%           gamma: Vector of gamma values returned from constraints
%                   [nc+nb x 1]
%               NOTE: Constraints MUST be in order of: all kinematic
%               constraints, then all euler parameterization constraints
%           phi_q: Jacobian of constraint functions [nc+nb x 7*nb]
% 
% Outputs:  q_dd: Accelerations in [r p]' format
%           lam:  Reaction Force Lagrange Multipliers
%               NOTE: Formatted similarly to gamma [lam lam_p]' [nc+nb x 1]
%                     lam corresponds to kinematic constraints [nc x 1]
%                     lam_p corresponds to euler parameter constraints
%                     [nb x 1]
%
% Written By: Dan Piombino
% 11/3/16

nb = length(M)/3;
nc = length(gamma)-nb;

z12 = zeros(3*nb,4*nb);
z14 = zeros(3*nb,nb);
z21 = zeros(4*nb,3*nb);
z33 = zeros(nc,nc);
z34 = zeros(nc,nb);
z41 = zeros(nb,3*nb);
z43 = zeros(nb,nc);
z44 = zeros(nb,nb);

phi_r = phi_q(1:nc,1:3*nb);
phi_p = phi_q(1:nc,3*nb+1:end);
P = phi_q(nc+1:end,3*nb+1:end);

LHS = [M z12 phi_r' z14;z21 Jp phi_p' P';phi_r phi_p z33 z34;z41 P z43 z44];
RHS = [F;tau;gamma(end);gamma(1:end-1)];
z = inv(LHS)*RHS;
[q_dd,lam] = decz(z,nb);
end