% simEngine3D_A8P1.m
% creates model for HW6 Problem 3
%
% Written by: Dan Piombino
% 11/3/16

clear;
clc;
tic;

%% Initial Conditions and Constants
%Environment
g = 9.81*[0 0 -1]'; %gravity in m/s^2

%System
nb = 2; %2 bodies
nc = 10; %10 non-ep constraints
h = 1e-4; %time step
T = 10; %length of simulation

%Body 1
L1 = 2; %half length of rod1
w1 = 0.05; %x-section width of rod1
Area1 = w1^2; %area of rod1
rho1 = 7800; %density of rod1
m1 = Area1*L1*rho1; %mass of rod1
M1 = eye(3)*m1; %Mass matrix of rod1
J1_bx = m1/12*(w1^2+w1^2);
J1_by = m1/12*(w1^2+(2*L1)^2);
J1_bz = J1_by;
J1_b = [J1_bx 0 0; 0 J1_by 0; 0 0 J1_bz];
Fm1 = m1*g; %weight of rod1

th1i = pi/2; %initial angle of rod1 (w/ respect to global -z)
r1_i = [0 L1 0]'; %position of point O1 at time 0
A1_i = [0 0 1;sin(th1i) cos(th1i) 0; -cos(th1i) sin(th1i) 0];
p1_i = iorient(A1_i); %euler parameters of rod1 initial rotation
q1_d_i = zeros(7,1);
q1_dd_i = zeros(7,1);

Fa1 = [0 0 0]';
nm1_b = [0 0 0]';
na1_b = [0 0 0]';
n1_b = na1_b+nm1_b;
F1 = Fa1+Fm1;

%Body 2
L2 = 1; %half length of rod2
w2 = 0.05; %x-section width of rod2
Area2 = w2^2; %area of rod2
rho2 = 7800; %density of rod2
m2 = Area2*L2*rho2; %mass of rod2
M2 = eye(3)*m2; %Mass matrix of rod2
J2_bx = m2/12*(w2^2+w2^2);
J2_by = m2/12*(w2^2+(2*L2)^2);
J2_bz = J2_by;
J2_b = [J2_bx 0 0; 0 J2_by 0; 0 0 J2_bz];
Fm2 = m2*g; %weight of rod2

th2i = 0; %initial angle of rod2 (w/ respect to global -z)
r2_i = [0 2*L1 -L2]'; %position of point O2 at time 0
A2_i = [0 0 1;0 1 0;-1 0 0];
p2_i = iorient(A2_i); %euler parameters of rod2 initial rotation
q2_d_i = zeros(7,1);
q2_dd_i = zeros(7,1);

Fa2 = [0 0 0]';
nm2_b = [0 0 0]';
na2_b = [0 0 0]';
n2_b = na2_b+nm2_b;
F2 = Fa2+Fm2;

%Combining Bodies
q(:,1) = [r1_i;r2_i;p1_i;p2_i;];
q_d(:,1) = [q1_d_i;q2_d_i];
q_dd(:,1) = [q1_dd_i;q2_dd_i]';
M = [M1 zeros(3,3);zeros(3,3) M2]; %combines mass matrices

%Reference Vectors
ax = [1 0 0]';
ay = [0 1 0]';
az = [0 0 1]';
sq1_0b = [0 0 0]';
sq1_1b = [-L1 0 0]';
sq2_1b = [L1 0 0]';
sq2_2b = [-L2 0 0]';
dis = [0 0 0]'; %reference distance for cd's
th = [pi/2 0 0]'; %reference angle for dp1's

%% Evaluate Constraints @ t=0 (NR to ensure consistent IC's)

nt = 1;
k = 1;
n = 1;
while n == 1
    %Kinematic Constraints
    [phi(1,k),nu(1,1),gamma(1,1),phi_r(1,:),phi_p(1,:)] = phi_cd(q(:,nt),q_d(:,nt),1,0,sq1_1b,sq1_0b,ax,dis);
    [phi(2,k),nu(2,1),gamma(2,1),phi_r(2,:),phi_p(2,:)] = phi_cd(q(:,nt),q_d(:,nt),1,0,sq1_1b,sq1_0b,ay,dis);
    [phi(3,k),nu(3,1),gamma(3,1),phi_r(3,:),phi_p(3,:)] = phi_cd(q(:,nt),q_d(:,nt),1,0,sq1_1b,sq1_0b,az,dis);
    [phi(4,k),nu(4,1),gamma(4,1),phi_r(4,:),phi_p(4,:)] = phi_cd(q(:,nt),q_d(:,nt),2,1,sq2_2b,sq2_1b,ax,dis);
    [phi(5,k),nu(5,1),gamma(5,1),phi_r(5,:),phi_p(5,:)] = phi_cd(q(:,nt),q_d(:,nt),2,1,sq2_2b,sq2_1b,ay,dis);
    [phi(6,k),nu(6,1),gamma(6,1),phi_r(6,:),phi_p(6,:)] = phi_cd(q(:,nt),q_d(:,nt),2,1,sq2_2b,sq2_1b,az,dis);
    [phi(7,k),nu(7,1),gamma(7,1),phi_r(7,:),phi_p(7,:)] = phi_dp1(q(:,nt),q_d(:,nt),1,0,az,az,th);
    [phi(8,k),nu(8,1),gamma(8,1),phi_r(8,:),phi_p(8,:)] = phi_dp1(q(:,nt),q_d(:,nt),1,0,az,ay,th);
    [phi(9,k),nu(9,1),gamma(9,1),phi_r(9,:),phi_p(9,:)] = phi_dp1(q(:,nt),q_d(:,nt),2,0,az,az,th);
    [phi(10,k),nu(10,1),gamma(10,1),phi_r(10,:),phi_p(10,:)] = phi_dp1(q(:,nt),q_d(:,nt),2,0,az,ay,th);
    [phi(11,k),nu(11,1),gamma(11,1),phi_r(11,:),phi_p(11,:)] = phi_cd(q(:,nt),q_d(:,nt),1,0,[0 0 0]',[1 L1 0]',ay,dis); %at t=0
    [phi(12,k),nu(12,1),gamma(12,1),phi_r(12,:),phi_p(12,:)] = phi_cd(q(:,nt),q_d(:,nt),2,0,[0 0 0]',[0 2*L1 -L2]',az,dis); %at t=0
    %EP Constraints
    for i = 1:nb
        [phi(nc+2+i,k),nu(nc+2+i,1),gamma(nc+2+i,1),phi_r(nc+2+i,:),phi_p(nc+2+i,:)] = phi_ep(q(:,nt),q_d(:,nt),i);
    end
    phi_q = [phi_r phi_p]; %Combined Jacobian
    
    dum1 = 1;
    dum2 = 0;
    for n2 = 1:size(phi,1);
        if abs(phi(n2,k)) >= 1e-8 %are any phi values still not zero?
            dum1 = 0; %if so, do not stop while loop
        end
        if k > 1
            if abs(phi(n2,k)-phi(n2,k-1)) >= 1e-10 %are any phi values still changing significantly?
                dum2 = 0; %if so, do not stop while loop
            end
        else
            dum2 = 0;
        end
    end
    if dum1 == 1 %if every value of phi is below threshold
        n = 0;
        iter(nt) = k;
    elseif dum2 == 1 %if every value of phi has converged to SOMETHING
        n = 0;
        itern(nt) = k;
    elseif k > 1000 %if iteration limit is reached
        error('The solver failed to converge')
    else
        q(:,nt) = q(:,nt)-inv(phi_q)*phi(:,k); %Newton-Raphson
        q_d(:,nt) = inv(phi_q)*nu; %Update q_dot
        q_dd(:,nt) = inv(phi_q)*gamma; %update q_doubledot
        k = k+1;
    end
end

clear phi phi_q phi_p phi_r nu gamma
%Calculate constraint info without additional starting conditions
[phi(1,1),nu(1,1),gamma(1,1),phi_r(1,:),phi_p(1,:)] = phi_cd(q(:,nt),q_d(:,nt),1,0,sq1_1b,sq1_0b,ax,dis);
[phi(2,1),nu(2,1),gamma(2,1),phi_r(2,:),phi_p(2,:)] = phi_cd(q(:,nt),q_d(:,nt),1,0,sq1_1b,sq1_0b,ay,dis);
[phi(3,1),nu(3,1),gamma(3,1),phi_r(3,:),phi_p(3,:)] = phi_cd(q(:,nt),q_d(:,nt),1,0,sq1_1b,sq1_0b,az,dis);
[phi(4,1),nu(4,1),gamma(4,1),phi_r(4,:),phi_p(4,:)] = phi_cd(q(:,nt),q_d(:,nt),2,1,sq2_2b,sq2_1b,ax,dis);
[phi(5,1),nu(5,1),gamma(5,1),phi_r(5,:),phi_p(5,:)] = phi_cd(q(:,nt),q_d(:,nt),2,1,sq2_2b,sq2_1b,ay,dis);
[phi(6,1),nu(6,1),gamma(6,1),phi_r(6,:),phi_p(6,:)] = phi_cd(q(:,nt),q_d(:,nt),2,1,sq2_2b,sq2_1b,az,dis);
[phi(7,1),nu(7,1),gamma(7,1),phi_r(7,:),phi_p(7,:)] = phi_dp1(q(:,nt),q_d(:,nt),1,0,az,az,th);
[phi(8,1),nu(8,1),gamma(8,1),phi_r(8,:),phi_p(8,:)] = phi_dp1(q(:,nt),q_d(:,nt),1,0,az,ay,th);
[phi(9,1),nu(9,1),gamma(9,1),phi_r(9,:),phi_p(9,:)] = phi_dp1(q(:,nt),q_d(:,nt),2,0,az,az,th);
[phi(10,1),nu(10,1),gamma(10,1),phi_r(10,:),phi_p(10,:)] = phi_dp1(q(:,nt),q_d(:,nt),2,0,az,ay,th);
%EP Constraints
for i = 1:nb
    [phi(nc+i,1),nu(nc+i,1),gamma(nc+i,1),phi_r(nc+i,:),phi_p(nc+i,:)] = phi_ep(q(:,nt),q_d(:,nt),i);
end
phi_q = [phi_r phi_p]; %Combined Jacobian

%% "ma=F" for initial (t = 0) q_dd values

%Expanded "M" matrix (Only Jp is unknown)
G1 = gmat(q(:,nt),1);
G1_d = gmat(q_d(:,nt),1);
G2 = gmat(q(:,nt),2);
G2_d = gmat(q_d(:,nt),2);
J1p = 4*G1'*J1_b*G1;
J2p = 4*G2'*J2_b*G2;
Jp = [J1p zeros(4,4);zeros(4,4) J2p]; %compiled J matrix

%Forces & Torques
F = [F1;F2];
tau1_b = 2*G1'*n1_b+8*G1_d'*J1_b*G1_d*q([3*nb+1:3*nb+4],nt);
tau2_b = 2*G2'*n2_b+8*G2_d'*J1_b*G2_d*q([4+3*nb+1:4+3*nb+4],nt);
tau_b = [tau1_b;tau2_b];

%Solve for q_dd and lamda
[q_dd(:,nt),lam(:,nt)] = maf(M,Jp,phi_q,F,tau_b,gamma);

%Reaction forces @ time t=0
R = -phi_q'*lam(:,nt);
Fr_O1(:,nt) = R(1:3); %Global reaction forces at O1
Fr_O2(:,nt) = R(4:6); %Global reaction forces at O2
A1 = orient(q(:,nt),1);
A2 = orient(q(:,nt),2);
nr_O1(:,nt) = A1*0.5*G1*R(7:10); %Global reaction torques at O1
nr_O2(:,nt) = A2*0.5*G2*R(11:14); %Global reaction torques at O2

%Time Response @ t=0
O1(:,nt) = q(1:3,nt);
O1_d(:,nt) = q_d(1:3,nt);
O1_dd(:,nt) = q_dd(1:3,nt);
O2(:,nt) = q(4:6,nt);
O2_d(:,nt) = q_d(4:6,nt);
O2_dd(:,nt) = q_dd(4:6,nt);

B1 = bmat(q(:,nt),1,sq1_1b);
B1_d = bmat(q_d(:,nt),1,sq1_1b);
Q(:,nt) = O1(:,nt)+A1*sq1_1b; %Location of end of rod (should be origin)
Q_d(:,nt) = O1_d(:,nt)+B1*q_d(7:10,nt);
Q_dd(:,nt) = O1_dd(:,nt)+B1_d*q_d(7:10,nt)+B1*q_dd(7:10,nt);

%% Dynamic Analysis
t(1) = 0;
for nt = 2:T/h+1
    t(nt) = (nt-1)*h;
    
    %Initial Guess for time nt
    q_dd(:,nt) = q_dd(:,nt-1);
    lam(:,nt) = lam(:,nt-1);
    
    n = 1;
    k = 1;
    clear res norm
    while n == 1
        
        %BDF to get q_d and q
%         clear r_dd p_dd phi phi_q phi_r phi_p psi nu gamma
%         if nt < 3 %not enough for BDF 2
            C = 1;
            beta = 1;
            q_d(:,nt) = C*q_d(:,nt-1)+beta*h*q_dd(:,nt);
            q(:,nt) = C*q(:,nt-1)+beta*h*q_d(:,nt);
%         else %BDF 2
%             C1 = 4/3; %1 point back
%             C2 = -1/3; %2 points back
%             beta = 2/3;
%             q_d(:,nt) = C1*q_d(:,nt-1)+C2*q_d(:,nt-2)+beta*h*q_dd(:,nt);
%             q(:,nt) = C1*q(:,nt-1)+C2*q(:,nt-2)+beta*h*q_d(:,nt);
%         end
        
        %% Calculate g(z)
        
        %Kinematic Constraints
        [phi(1,k),nu(1,1),gamma(1,1),phi_r(1,:),phi_p(1,:)] = phi_cd(q(:,nt),q_d(:,nt),1,0,sq1_1b,sq1_0b,ax,dis);
        [phi(2,k),nu(2,1),gamma(2,1),phi_r(2,:),phi_p(2,:)] = phi_cd(q(:,nt),q_d(:,nt),1,0,sq1_1b,sq1_0b,ay,dis);
        [phi(3,k),nu(3,1),gamma(3,1),phi_r(3,:),phi_p(3,:)] = phi_cd(q(:,nt),q_d(:,nt),1,0,sq1_1b,sq1_0b,az,dis);
        [phi(4,k),nu(4,1),gamma(4,1),phi_r(4,:),phi_p(4,:)] = phi_cd(q(:,nt),q_d(:,nt),2,1,sq2_2b,sq2_1b,ax,dis);
        [phi(5,k),nu(5,1),gamma(5,1),phi_r(5,:),phi_p(5,:)] = phi_cd(q(:,nt),q_d(:,nt),2,1,sq2_2b,sq2_1b,ay,dis);
        [phi(6,k),nu(6,1),gamma(6,1),phi_r(6,:),phi_p(6,:)] = phi_cd(q(:,nt),q_d(:,nt),2,1,sq2_2b,sq2_1b,az,dis);
        [phi(7,k),nu(7,1),gamma(7,1),phi_r(7,:),phi_p(7,:)] = phi_dp1(q(:,nt),q_d(:,nt),1,0,az,az,th);
        [phi(8,k),nu(8,1),gamma(8,1),phi_r(8,:),phi_p(8,:)] = phi_dp1(q(:,nt),q_d(:,nt),1,0,az,ay,th);
        [phi(9,k),nu(9,1),gamma(9,1),phi_r(9,:),phi_p(9,:)] = phi_dp1(q(:,nt),q_d(:,nt),2,0,az,az,th);
        [phi(10,k),nu(10,1),gamma(10,1),phi_r(10,:),phi_p(10,:)] = phi_dp1(q(:,nt),q_d(:,nt),2,0,az,ay,th);
        %EP Constraints
        for i = 1:nb
            [phi(nc+i,k),nu(nc+i,1),gamma(nc+i,1),phi_r(nc+i,:),phi_p(nc+i,:)] = phi_ep(q(:,nt),q_d(:,nt),i);
        end
        phi_q = [phi_r phi_p]; %Combined Jacobian
        
        %Expanded "M" matrix (Only Jp is unknown)
        G1 = gmat(q(:,nt),1);
        G1_d = gmat(q_d(:,nt),1);
        G2 = gmat(q(:,nt),2);
        G2_d = gmat(q_d(:,nt),2);
        J1p = 4*G1'*J1_b*G1;
        J2p = 4*G2'*J2_b*G2;
        Jp = [J1p zeros(4,4);zeros(4,4) J2p]; %compiled J matrix
        
        %Forces & Torques
        F = [F1;F2];
        tau1_b = 2*G1'*n1_b+8*G1_d'*J1_b*G1_d*q([3*nb+1:3*nb+4],nt);
        tau2_b = 2*G2'*n2_b+8*G2_d'*J2_b*G2_d*q([4+3*nb+1:4+3*nb+4],nt);
        tau_b = [tau1_b;tau2_b];
        
        %Calculate residual (g(z))
        res(:,k) = bdfg(M,Jp,q_dd(:,nt),F,tau_b,h,beta,phi(:,k),phi_q,lam(:,nt)); %calculates residual
        
        %% Convergeance Logic
        
        dum1 = 1;
        dum2 = 1;
        for n2 = 1:size(res,1);
            if abs(res(n2,k)) >= 1e-3 %are any phi values still not zero?
                dum1 = 0; %if so, do not stop while loop
            end
            if k > 1
                if abs(res(n2,k)-res(n2,k-1)) >= 1e-6 %are any phi values still changing significantly?
                    dum2 = 0; %if so, do not stop while loop
                end
            else
                dum2 = 0;
            end
        end
        if dum1 == 1 %if every value of res is below threshold
            n = 0;
            iter(nt) = k;
        elseif dum2 == 1 %if every value of res has converged to SOMETHING
            n = 0;
            iter(nt) = k;
        elseif k > 10000 %if iteration limit is reached
            error('The solver failed to converge')
        else %if not converged
            %Compute psi
            psi = qnpsi(M,Jp,phi_q);
            z = [q_dd(:,nt);lam(:,nt)];
            z = z-inv(psi)*res(:,k);
            [q_dd(:,nt),lam(:,nt)] = decz(z,nb);
            k = k+1;
            clear r_dd p_dd
        end
    end %end of while loop
    %% Desired Quantities @ step nt
    %Velocity Constraints
    vc(:,nt) = phi_q*q_d(:,nt)-nu;
    nvc(nt) = norm(vc([4 5 6 9 10],nt));
    
    %Reaction forces & Torques
    R = -phi_q'*lam(:,nt);
    Fr_O1(:,nt) = R(1:3); %Global reaction forces at O1
    Fr_O2(:,nt) = R(4:6); %Global reaction forces at O2
    A1 = orient(q(:,nt),1);
    A2 = orient(q(:,nt),2);
    nr_O1(:,nt) = A1*0.5*G1*R(7:10); %Global reaction torques at O1
    nr_O2(:,nt) = A2*0.5*G2*R(11:14); %Global reaction torques at O2
    
    %Time Response @ t=0
    O1(:,nt) = q(1:3,nt);
    O1_d(:,nt) = q_d(1:3,nt);
    O1_dd(:,nt) = q_dd(1:3,nt);
    O2(:,nt) = q(4:6,nt);
    O2_d(:,nt) = q_d(4:6,nt);
    O2_dd(:,nt) = q_dd(4:6,nt);
    
    B1 = bmat(q(:,nt),1,sq1_1b);
    B1_d = bmat(q_d(:,nt),1,sq1_1b);
    Q(:,nt) = O1(:,nt)+A1*sq1_1b; %Location of end of rod (should be origin)
    Q_d(:,nt) = O1_d(:,nt)+B1*q_d(7:10,nt);
    Q_dd(:,nt) = O1_dd(:,nt)+B1_d*q_d(7:10,nt)+B1*q_dd(7:10,nt);
end

trun = toc

%% Plot Results
% Point O1
figure
hold on
plot(t,O1(1,:))
plot(t,O1(2,:))
plot(t,O1(3,:))
title('Position of Point O1-prime')
ylabel('Position [m]')
xlabel('Time [s]')
legend('Global X','Global Y','Global Z')

figure
hold on
plot(t,O1_d(1,:))
plot(t,O1_d(2,:))
plot(t,O1_d(3,:))
title('Velocity of Point O1-prime')
ylabel('Velocity [m/s]')
xlabel('Time [s]')
legend('Global X','Global Y','Global Z')

figure
hold on
plot(t,O1_dd(1,:))
plot(t,O1_dd(2,:))
plot(t,O1_dd(3,:))
title('Acceleration of Point O1-prime')
ylabel('Acceleration [m/s^2]')
xlabel('Time [s]')
legend('Global X','Global Y','Global Z')

%Point O2
figure
hold on
plot(t,O2(1,:))
plot(t,O2(2,:))
plot(t,O2(3,:))
title('Position of Point O2-prime')
ylabel('Position [m]')
xlabel('Time [s]')
legend('Global X','Global Y','Global Z')

figure
hold on
plot(t,O2_d(1,:))
plot(t,O2_d(2,:))
plot(t,O2_d(3,:))
title('Velocity of Point O2-prime')
ylabel('Velocity [m/s]')
xlabel('Time [s]')
legend('Global X','Global Y','Global Z')

figure
hold on
plot(t,O2_dd(1,:))
plot(t,O2_dd(2,:))
plot(t,O2_dd(3,:))
title('Acceleration of Point O2-prime')
ylabel('Acceleration [m/s^2]')
xlabel('Time [s]')
legend('Global X','Global Y','Global Z')

% %Point Q
% figure
% hold on
% plot(t,Q(1,:))
% plot(t,Q(2,:))
% plot(t,Q(3,:))
% title('Position of Point Q')
% ylabel('Position [m]')
% xlabel('Time [s]')
% legend('Global X','Global Y','Global Z')
% 
% figure
% hold on
% plot(t,Q_d(1,:))
% plot(t,Q_d(2,:))
% plot(t,Q_d(3,:))
% title('Velocity of Point Q')
% ylabel('Velocity [m/s]')
% xlabel('Time [s]')
% legend('Global X','Global Y','Global Z')
% 
% figure
% hold on
% plot(t,Q_dd(1,:))
% plot(t,Q_dd(2,:))
% plot(t,Q_dd(3,:))
% title('Acceleration of Point Q')
% ylabel('Acceleration [m/s^2]')
% xlabel('Time [s]')
% legend('Global X','Global Y','Global Z')

figure
hold on
plot(t,nvc)
title('Revolute Joint 2 Velocity Constraint Violation')
ylabel('Norm of Violation')
xlabel('Time [s]')
