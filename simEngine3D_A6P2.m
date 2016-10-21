% simEngine3D_A6P2.m
% creates model for HW6 Problem 2
%
% Written by: Dan Piombino
% 10/20/16

clear;
clc;

%% Create Model
t = 0;
theta(1,1) = pi/4*cos(2*t);
theta(2,1) = -pi/2*sin(2*t);
theta(3,1) = -pi*cos(2*t);
thg = 40*pi/180; %theta "guess" @ 40 deg instead of actual 45
r0 = [0 2*sin(thg) -2*cos(thg)]'; %initial r configuration guess
A0 = [0 0 1;sin(thg) cos(thg) 0; -cos(thg) sin(thg) 0];
p0 = iorient(A0); %initial orientation guess
q = [r0;p0];
q_d = zeros(7,1); %initial velocity guess
L = 2;
i = 1;
j = 0;
sq_ib = [-L,0,0]';
so_jb = [0 0 0]';
ax_jb = [1 0 0]';
ax_ib = [1 0 0]';
az_ib = [0 0 1]';
az_jb = [0 0 1]';
ay_jb = [0 1 0]';
c1 = [1 0 0]';
c2 = [0 1 0]';
c3 = [0 0 1]';
dis = zeros(3,1);
th2 = [pi/2 0 0]';

n = 1;
k = 1;
while n == 1
    [phi(1,k),~,~,phi_r(1,:),phi_p(1,:)] = phi_dp1(q,q_d,i,j,az_ib,ay_jb,th2); 
    [phi(2,k),~,~,phi_r(2,:),phi_p(2,:)] = phi_cd(q,q_d,i,j,sq_ib,so_jb,c1,dis); %fixed in x
    [phi(3,k),~,~,phi_r(3,:),phi_p(3,:)] = phi_cd(q,q_d,i,j,sq_ib,so_jb,c2,dis); %fixed in y
    [phi(4,k),~,~,phi_r(4,:),phi_p(4,:)] = phi_cd(q,q_d,i,j,sq_ib,so_jb,c3,dis); %fixed in z
    [phi(5,k),~,~,phi_r(5,:),phi_p(5,:)] = phi_dp1(q,q_d,i,j,az_ib,az_jb,th2);
    [phi(6,k),~,~,phi_r(6,:),phi_p(6,:)] = phi_dp1(q,q_d,i,j,ax_ib,-az_jb,theta);
    [phi(7,k),~,~,phi_r(7,:),phi_p(7,:)] = phi_ep(q,q_d,i);
    phi_q = [phi_r phi_p];
    
    dum1 = 1;
    dum2 = 1;
    for n2 = 1:size(phi,1);
        if abs(phi(n2,k)) >= 0.001
            dum1 = 0;
        end
        if k > 1
            if abs(phi(n2,k)-phi(n2,k-1)) >= 1e-4
                dum2 = 0;
            end
        else
            dum2 = 0;
        end
    end
    if dum1 == 1 %if every value of phi is below 0.001
        n = 0;
        disp('The solution converged')
    elseif dum2 == 1
        n = 0;
        disp('The value of phi was not changing with further iterations')
    elseif k > 1000
        disp('The solver failed to converge')
        n = 0;
    else
        q = q-(phi_q^-1)*phi(:,k);
        k = k+1;
    end
end

[phi(1,k),nu(1,1),gamma(1,1),phi_r(1,:),phi_p(1,:)] = phi_dp1(q,q_d,i,j,az_ib,ax_jb,[0 0 0]'); 
[phi(2,k),nu(2,1),gamma(2,1),phi_r(2,:),phi_p(2,:)] = phi_cd(q,q_d,i,j,sq_ib,so_jb,c1,dis); %fixed in x
[phi(3,k),nu(3,1),gamma(3,1),phi_r(3,:),phi_p(3,:)] = phi_cd(q,q_d,i,j,sq_ib,so_jb,c2,dis); %fixed in y
[phi(4,k),nu(4,1),gamma(4,1),phi_r(4,:),phi_p(4,:)] = phi_cd(q,q_d,i,j,sq_ib,so_jb,c3,dis); %fixed in z
[phi(5,k),nu(5,1),gamma(5,1),phi_r(5,:),phi_p(5,:)] = phi_dp1(q,q_d,i,j,az_ib,az_jb,th2);
[phi(6,k),nu(6,1),gamma(6,1),phi_r(6,:),phi_p(6,:)] = phi_dp1(q,q_d,i,j,ax_ib,-az_jb,theta);
[phi(7,k),nu(7,1),gamma(7,1),phi_r(7,:),phi_p(7,:)] = phi_ep(q,q_d,i);
phi_q = [phi_r phi_p];

q_d = (phi_q^-1)*nu;
q_dd = (phi_q^-1)*gamma;

A = orient(q,i);
Q = q(1:3)+A*sq_ib

phif = phi(6,k)
phi_qf = phi_q(6,:)
nuf = nu(6)
gammaf = gamma(6)