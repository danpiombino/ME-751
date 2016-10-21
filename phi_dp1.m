function [phi,nu,gamma,phi_r,phi_p] = phi_dp1(q,q_d,i,j,ai_b,aj_b,theta)
% Calculates parameters associated with DP1 constraint between coordinate
% systems i and j
% 
% Inputs:   q = [7*nb x 1] vector of positions and euler parameters
%               in the form of [r1 r2 ... rnb p1 p2 ... pnb]'
%           i = Index of 1st coordinate system of interest (0 = ground)
%           j = Index of 2nd coordinate system of interest (0 = ground)
%           q_d = [7*nb x 1] vector of velocities in r and p directions
%           ai_b = DP1 reference vector from coordinate system i (local)
%           aj_b = DP1 reference vector from coordinate system j (local)
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
% 10/13/16

%% Setup

nb = length(q)/7;

% Coordinate System i
Ai = orient(q,i);
ai = Ai*ai_b;
Bi = bmat(q,i,ai_b);
Bi_d = bmat(q_d,i,ai_b);

if i == 0
    pi = [1 0 0 0]'; %by definition
    pi_d = zeros(4,1); %by definition
else
    pi = q([4*(i-1)+3*nb+1:4*(i-1)+3*nb+4]);
    pi_d = q_d([4*(i-1)+3*nb+1:4*(i-1)+3*nb+4]);
end

% Coordinate System j
Aj = orient(q,j);
aj = Aj*aj_b;
Bj = bmat(q,j,aj_b);
Bj_d = bmat(q_d,j,aj_b);

if j == 0
    pj = [1 0 0 0]'; %by definition
    pj_d = zeros(4,1); %by definition
else
    pj = q([4*(j-1)+3*nb+1:4*(j-1)+3*nb+4]);
    pj_d = q_d([4*(j-1)+3*nb+1:4*(j-1)+3*nb+4]);
end

% Theta 
f = norm(ai)*norm(aj)*cos(theta(1));
f_d = -norm(ai)*norm(aj)*sin(theta(1))*theta(2);     
f_dd = -norm(ai)*norm(aj)*(cos(theta(1))*theta(2)^2+sin(theta(1))*theta(3));

%% Evaluate Quantities

phi = ai_b'*Ai'*Aj*aj_b-f;
nu = f_d;
gamma = -ai'*Bj_d*pj_d-aj'*Bi_d*pi_d-2*(Bi*pi_d)'*(Bj*pj_d)+f_dd;

phi_r = zeros(1,3*nb);  %by definition for DP1
phi_p = zeros(1,4*nb);
phi_pi = aj'*Bi; 
phi_pj = ai'*Bj;
if i ~= 0
    phi_p(4*(i-1)+1:4*(i-1)+4) = phi_pi;
end
if j ~= 0
    phi_p(4*(j-1)+1:4*(j-1)+4) = phi_pj;
end
phi_q = [phi_r phi_p];
end