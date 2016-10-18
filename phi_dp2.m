function [phi,nu,gamma,phi_r,phi_p] = phi_dp2(q,q_d,i,j,ai_b,sp_ib,sq_jb,theta)
% Calculates parameters associated with DP1 constraint between coordinate
% systems i and j
% 
% Inputs:   q = [7*nb x 1] vector of positions and euler parameters
%               in the form of [r1 r2 ... rnb p1 p2 ... pnb]'
%           i = Index of 1st coordinate system of interest (0 = ground)
%           j = Index of 2nd coordinate system of interest (0 = ground)
%           q_d = [7*nb x 1] vector of velocities in r and p directions
%           ai_b = DP1 reference vector from coordinate system i (local)
%           sp_ib = vector in local coordinate system i to reference point
%                   p: [3x1 column vector]
%           sq_ib = vector in local coordinate system j to reference point
%                   q: [3x1 column vector]
%           theta = Desired angle between reference vectors, theta, and 
%                   derivatives over time. Input should take the form:
%                   [theta(t) theta_dot(t) theta_doubledot(t)]'
%                   UNITS: [RAD RAD/S RAD/S/S]'
%
% Outputs:  phi = value of constraint evaluated at current conditions
%           nu = coefficeints on right side of velocity equation
%           gamma = coefficeints on right side of acceleration equation
%           phi_r = partial derivative of phi with respect to r
%           phi_p = partial derivative of phi with respect to p
%
% Written by: Dan Piombino
% 10/17/16

%% Setup

nb = length(q)/7;

% Coordinate System i
Ai = orient(q,i);
ai = Ai*ai_b;
Bia = bmat(q,i,ai_b);
Bia_d = bmat(q_d,i,ai_b);
Bip = bmat(q,i,sp_ib);
Bip_d = bmat(q_d,i,sp_ib);

if i == 0
    ri = [0 0 0]'; %by definition
    ri_d = [0 0 0]'; %by definition
    pi = [1 0 0 0]'; %by definition
    pi_d = zeros(4,1); %by definition
else
    ri = q([3*(i-1)+1:3*(i-1)+3]);
    ri_d = q_d([3*(i-1)+1:3*(i-1)+3]);
    pi = q([4*(i-1)+3*nb+1:4*(i-1)+3*nb+4]);
    pi_d = q_d([4*(i-1)+3*nb+1:4*(i-1)+3*nb+4]);
end

% Coordinate System j
Aj = orient(q,j);
Bjq = bmat(q,j,sq_jb);
Bjq_d = bmat(q_d,j,sq_jb);

if j == 0
    rj = [0 0 0]'; %by definition
    rj_d = [0 0 0]'; %by definition
    pj = [1 0 0 0]'; %by definition
    pj_d = zeros(4,1); %by definition
else
    rj = q([3*(i-1)+1:3*(i-1)+3]);
    rj_d = q_d([3*(i-1)+1:3*(i-1)+3]);
    pj = q([4*(j-1)+3*nb+1:4*(j-1)+3*nb+4]);
    pj_d = q_d([4*(j-1)+3*nb+1:4*(j-1)+3*nb+4]);
end

% Theta 
f = norm(ai_b)*norm(aj_b)*cos(theta(1));
f_d = -norm(ai_b)*norm(aj_b)*sin(theta(1))*theta(2);     
f_dd = -norm(ai_b)*norm(aj_b)*(cos(theta(1))*theta(2)^2+sin(theta(1))*theta(3));

%% Evaluate Quantities

phi = ai'*(rj+Aj*sq_jb-ri-Ai*sp_ib)-f;
nu = -ai'*(rj_d+Bjq*pj-ri_d-Bip*pi_d)-(rj+Aj*sq_jb-ri-Ai*sp_ib)'*Bia*pi_d+f_d; %nu = -phi_t
gamma = -ai'*Bjq_d*pj_d+ai'*Bip_d*pi_d-(rj+Aj*sq_jb-ri-Ai*sp_ib)'*Bia_d*pi_d-2*(Bia*pi_d)'*(rj_d+Bjq*pj-ri_d-Bip*pi_d)+f_dd;

phi_r = zeros(1,3*nb);
phi_ri = -ai';
phi_rj = ai';
phi_r(3*(i-1)+1:3*(i-1)+3) = phi_ri;
phi_r(3*(j-1)+1:3*(j-1)+3) = phi_rj;
phi_p = zeros(1,4*nb);
phi_pi = -ai'*Bip-sp_ib'*Ai'*Bia;
phi_pj = ai'*Bjq;
phi_p(4*(i-1)+1:4*(i-1)+4) = phi_pi;
phi_p(4*(j-1)+1:4*(j-1)+4) = phi_pj;
phi_q = [phi_r phi_p];
end