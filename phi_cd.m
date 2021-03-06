function [phi,nu,gamma,phi_r,phi_p,rf_q] = phi_cd(q,q_d,i,j,sp_ib,sq_jb,c,dis,lam)
% Calculates parameters associated with DP1 constraint between coordinate
% systems i and j
% 
% Inputs:   q = [7*nb x 1] vector of positions and euler parameters
%               in the form of [r1 r2 ... rnb p1 p2 ... pnb]'
%           i = Index of 1st coordinate system of interest (0 = ground)
%           j = Index of 2nd coordinate system of interest (0 = ground)
%           q_d = [7*nb x 1] vector of velocities in r and p directions
%           sp_ib = vector in local coordinate system i to reference point
%                   p: [3x1 column vector]
%           sq_ib = vector in local coordinate system j to reference point
%                   q: [3x1 column vector]
%           c = vector in global coordinate system along which to enforce
%               the constraint: [x y z]'
%           dis = column vector of prescribed displacement function and
%                 derivatives in the form of: [dis,dis_d,dis_dd]'
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
Bi = bmat(q,i,sp_ib);
Bi_d = bmat(q_d,i,sp_ib);

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
Bj = bmat(q,j,sq_jb);
Bj_d = bmat(q_d,j,sq_jb);

if j == 0
    rj = [0 0 0]'; %by definition
    rj_d = [0 0 0]'; %by definition
    pj = [1 0 0 0]'; %by definition
    pj_d = zeros(4,1); %by definition
else
    rj = q([3*(j-1)+1:3*(j-1)+3]);
    rj_d = q_d([3*(j-1)+1:3*(j-1)+3]);
    pj = q([4*(j-1)+3*nb+1:4*(j-1)+3*nb+4]);
    pj_d = q_d([4*(j-1)+3*nb+1:4*(j-1)+3*nb+4]);
end

%% Evaluate Quantities
phi = c'*(rj+Aj*sq_jb-ri-Ai*sp_ib)-dis(1);
% nu = -c'*(rj_d+Bj*pj_d-ri_d-Bi*pi_d)+dis(2);
nu = dis(2);
gamma = c'*(Bi_d*pi_d-Bj_d*pj_d)+dis(3);

phi_r = zeros(1,3*nb);
phi_ri = -c';
phi_rj = c';
phi_p = zeros(1,4*nb);
phi_pi = -c'*Bi;
phi_pj = c'*Bj;
if i ~= 0
    phi_r(3*(i-1)+1:3*(i-1)+3) = phi_ri;
    phi_p(4*(i-1)+1:4*(i-1)+4) = phi_pi;
end
if j ~= 0
    phi_r(3*(j-1)+1:3*(j-1)+3) = phi_rj;
    phi_p(4*(j-1)+1:4*(j-1)+4) = phi_pj;
end
phi_q = [phi_r phi_p];

if ~isempty(lam)
    rf_q = zeros(7*nb,7*nb);
    Ki = kmat(sp_ib,c);
    Kj = kmat(sq_jb,c);
    
    if i ~= 0
        rf_q(3*nb+4*(i-1)+1:3*nb+4*(i-1)+4,3*nb+4*(i-1)+1:3*nb+4*(i-1)+4) = -Ki;
    end
    if j ~= 0
        rf_q(3*nb+4*(j-1)+1:3*nb+4*(j-1)+4,3*nb+4*(j-1)+1:3*nb+4*(j-1)+4) = Kj;
    end
    rf_q = lam*rf_q;
else
    rf_q = [];
end
end