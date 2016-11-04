function res = bdfg(M,Jp,q_dd,F,tau,h,beta,phi,phi_q,lam)
% function res = bdfg(M,Jp,q_dd,F,tau,h,beta,phi,phi_q)
% Calculates the residual of g(z) at the given iteration
% NOTE: Constraints MUST be in order of: all kinematic
%       constraints, then all euler parameterization constraints
% 
% Inputs:   M: Mass matrix [3*nb x 3*nb]
%           Jp: Euler Parameter Inertia Matrix [4*nb x 4*nb]
%           q_dd: Accelerations in [r1 p1 r2 p2 ...]' format at desired
%           time
%           F: External Forces in LRF [3*nb x 1]
%           tau: External Torques in LRF [4*nb x 1]
%           h: BDF Method time step
%           phi: Evaluated constraint functions [nc+nb x 1]
%               NOTE: Constraints MUST be in order of: all kinematic
%               constraints, then all euler parameterization constraints
%           phi_q: Jacobian of constraint functions [nc+nb x 7*nb]
%           lam:  Reaction Force Lagrange Multipliers
%               NOTE: Formatted similarly to gamma [lam lam_p]' [nc+nb x 1]
%                     lam corresponds to kinematic constraints [nc x 1]
%                     lam_p corresponds to euler parameter constraints
%                     [nb x 1]
% 
% Outputs:  res: residual of g(z) function for bdf [8*nb+nc]
%
% Written By: Dan Piombino
% 11/3/16

nb = length(q_dd)/7;
nc = length(phi)-nb;
[r_dd,p_dd] = q2rp(q_dd);
phi_r = phi_q(1:nc,1:3*nb);
phi_p = phi_q(1:nc,3*nb+1:end);
P = phi_q(nc+1:end,3*nb+1:end);
lambda = lam(1:nc);
lamp = lam(nc+1:end);
a = (1/(h^2)*(beta^2));

res = [M*r_dd+phi_r'*lambda-F;Jp*p_dd+phi_p'*lambda+P'*lamp-tau;a*phi];
end




