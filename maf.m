function [q_dd,lam] = maf(M,Jp,phi_q,F,tau,gamma)
% function [q_dd,lam] = maf(M,Jp,phi_q,F,tau,gamma)
% calculates "a = inv(m)*F"
% Used for euqilibrium or inverse dynamics problems
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

LHS = zeros(8*nb+nc,8*nb+nc);
LHS(1:3*nb,1:3*nb) = M;
LHS(3*nb+1:7*nb,3*nb+1:7*nb) = Jp;
LHS(7*nb+1:end,1:7*nb) = phi_q;
LHS(1:7*nb,7*nb+1:end) = phi_q';
RHS = [F;tau;gamma];
z = inv(LHS)*RHS;
[q_dd,lam] = decz(z,nb);
end