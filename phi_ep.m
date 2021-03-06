function [phi,nu,gamma,phi_r,phi_p,rf_p] = phi_ep(q,q_d,i,lam)
% Calculates parameters associated with euler parameter constraint
% 
% Inputs:   q = [7*nb x 1] vector of positions and euler parameters
%               in the form of [r1 r2 ... rnb p1 p2 ... pnb]'
%           i = Index of coordinate system of interest
%           q_d = [7*nb x 1] vector of velocities in r and p directions
%           lam = lagrange multipliers [7nb+nc x 1]
%
% Outputs:  phi = value of constraint evaluated at current conditions
%           nu = coefficeints on right side of velocity equation
%           gamma = coefficeints on right side of acceleration equation
%           phi_r = partial derivative of phi with respect to r
%           phi_p = partial derivative of phi with respect to p4
%           rf_q = partial derivative of [phi_q*lam] with respect to q
%
% Written by: Dan Piombino
% 10/20/16

%% Setup

nb = length(q)/7;
p = q([4*(i-1)+3*nb+1:4*(i-1)+3*nb+4]);
p_d = q_d([4*(i-1)+3*nb+1:4*(i-1)+3*nb+4]);


%% Evaluate Quantities

phi = 0.5*(p'*p-1);
nu = 0; %nu = -phi_t
gamma = -p_d'*p_d;

phi_r = zeros(1,3*nb);  %by definition
phi_p = zeros(1,4*nb);
phi_pi = p';
if i ~= 0
    phi_p(4*(i-1)+1:4*(i-1)+4) = phi_pi;
end
phi_q = [phi_r phi_p];

if ~isempty(lam)
    rf_p = zeros(4*nb,4*nb);
    rf_p(4*(i-1)+1:4*(i-1)+4,i) = lam;
else
    rf_p = [];
end
end