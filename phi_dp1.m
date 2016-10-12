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
% 10/12/16

%% Setup

nb = length(q)/7;

% Coordinate System i
Ai = orient(q,i);
Ai_d = orientd(q,q_d,i);
wi_b = pdtowb(q,q_d,i);
Bi = bmat(q,i,ai_b);
ai_b_til = skew3(ai_b);
wi_b_til = skew3(wi_b);

% Coordinate System j
Aj = orient(q,j);
Aj_d = orientd(q,q_d,j);
wj_b = pdtowb(q,q_d,j);
Bj = bmat(q,j,aj_b);
aj_b_til = skew3(aj_b);
wj_b_til = skew3(wj_b);

% Theta 
f = norm(ai_b)*norm(aj_b)*cos(theta(1));
f_d = -norm(ai_b)*norm(aj_b)*sin(theta(2));     %unsure
f_dd = -norm(ai_b)*norm(aj_b)*cos(theta(3));    %unsure


%% Evaluate Quantities

phi = ai_b'*Ai'*Aj*aj_b-f;
nu = -(ai_b'*Ai'*aj_b'*Aj_d'+aj_b'*Aj'*ai_b'*Ai_d')+f_d;
gamma = -aj_b'*(Aj'*Ai*wi_b_til*wi_b_til+wj_b_til*wj_b_til*Aj'*Ai)*ai_b+2*wj_b'*aj_b_til*Aj*Ai*ai_b_til*wi_b+f_dd;
phi_r = zeros(1,3*nb);  %by definition for DP1
phi_p = zeros(1,4*nb);
phi_pi = aj_b'*Aj'*Bi;
phi_pj = ai_b'*Ai'*Bj;
phi_p(4*(i-1)+1:4*(i-1)+4) = phi_pi;
phi_p(4*(j-1)+1:4*(j-1)+4) = phi_pj;

end