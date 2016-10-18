% simEngine3D_A6P1.m
% runs phi_dp2 and phi_d functions for HW6 p-1
%
% Notes:    _d: refers to derivative
%           _dd: refers to double derivative
%           _ib: refers to local coordinate system i (read: i-bar)
%
% Written by: Dan Piombino
% 10/17/16

%% DP2

q = zeros(14,1); %[r1 r2 .... p1 p2 ...]'
q_d = zeros(14,1); %[r1_dot r2_dot ... p1_dot p2_dot ...]'
sp_ib1 = [1 1 1]'; %[x y z]' vector of point p in local reference frame i
sq_jb1 = [1 1 1]'; %[x y z]' vector of point q in local reference frame j
ai_b = [.5774 .5774 .5774]'; %[x y z]' vector of reference vector a on LRF i (does not need to be unit length)
theta = [0 0 0]'; %[theta(t) theta_dot(t) theta_doubledot(t)]' IN RADIANS
i1 = 1; %index of body i (0 = ground)
j1 = 0; %index of body j (0 = ground)

[phi1,nu1,gamma1,phi_r1,phi_p1] = phi_dp2(q,q_d,i1,j1,ai_b,sp_ib1,sq_jb1,theta);

%% Setup: D

sp_ib2 = [1 1 1]'; %[x y z]' vector of point p in local reference frame i
sq_jb2 = [1 1 1]'; %[x y z]' vector of point q in local reference frame j
dis = [0 0 0]'; %[distance velocity acceleration]' prescribed by constraint
i2 = 1; %index of body i (0 = ground)
j2 = 0; %index of body j (0 = ground)

[phi2,nu2,gamma2,phi_r2,phi_p2] = phi_d(q,q_d,i2,j2,sp_ib2,sq_jb2,dis);