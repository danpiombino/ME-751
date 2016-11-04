% simEngine3D_A8P1.m
% creates model for HW6 Problem 3
%
% Written by: Dan Piombino
% 11/3/16

clear;
clc;
tic;

%% Initial Conditions and Constants
thg = 40*pi/180; %theta "guess" @ 40 deg instead of actual 45
r0 = [0 2*sin(thg) -2*cos(thg)]'; %initial r configuration guess
A0 = [0 0 1;sin(thg) cos(thg) 0; -cos(thg) sin(thg) 0];
p0 = iorient(A0); %initial orientation guess
q(:,1) = [r0;p0]; %initial variables
q_d(:,1) = zeros(7,1); %initial velocity guess
q_dd(:,1) = zeros(7,1);
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
ar_jb = [0 -1 -1]';
c1 = [1 0 0]';
c2 = [0 1 0]';
c3 = [0 0 1]';
dis = zeros(3,1);
th2 = [pi/2 0 0]';
h = 1e-3 %time step
T = 10; %length of simulation
nb = 1;

%% Mass and other properties

rho = 7800;
w = .05; %x-section width
A = w^2; %area
m = A*L*rho;
M = eye(3)*m; %since there's only 1 body
Fm = m*9.81*[0 0 -1]'; %Force due to gravity on body i in global reference frame

J1_bx = m/12*(w^2+w^2);
J1_by = m/12*(w^2+(2*L)^2);
J1_bz = J1_by;
J1_b = [J1_bx 0 0; 0 J1_by 0; 0 0 J1_bz];

%% Newton Raphson on Initial Condition Guess

n = 1;
k = 1;
nt = 1; %time t=0
t(nt) = (nt-1)*h;
theta(1,nt) = pi/4*cos(2*t(nt))-pi/2;
theta(2,nt) = -pi/2*sin(2*t(nt));
theta(3,nt) = -pi*cos(2*t(nt));
    while n == 1
        %Kinematic Constraints
        [phi(1,k),nu(1,1),gamma(1,1),phi_r(1,:),phi_p(1,:)] = phi_dp1(q(:,nt),q_d(:,nt),i,j,az_ib,ay_jb,th2);
        [phi(2,k),nu(2,1),gamma(2,1),phi_r(2,:),phi_p(2,:)] = phi_cd(q(:,nt),q_d(:,nt),i,j,sq_ib,so_jb,c1,dis); %fixed in x
        [phi(3,k),nu(3,1),gamma(3,1),phi_r(3,:),phi_p(3,:)] = phi_cd(q(:,nt),q_d(:,nt),i,j,sq_ib,so_jb,c2,dis); %fixed in y
        [phi(4,k),nu(4,1),gamma(4,1),phi_r(4,:),phi_p(4,:)] = phi_cd(q(:,nt),q_d(:,nt),i,j,sq_ib,so_jb,c3,dis); %fixed in z
        [phi(5,k),nu(5,1),gamma(5,1),phi_r(5,:),phi_p(5,:)] = phi_dp1(q(:,nt),q_d(:,nt),i,j,az_ib,az_jb,th2);
        %Driving Constraints
        if abs(theta(1,nt))*180/pi >= 10 %Change reference angle to avoid singularities
            [phi(6,k),nu(6,1),gamma(6,1),phi_r(6,:),phi_p(6,:)] = phi_dp1(q(:,nt),q_d(:,nt),i,j,ax_ib,ay_jb,theta(:,nt)); %-az_jb
        else
            thetadev = theta(:,nt);
            thetadev(1) = theta(1,nt)+pi/4;
            [phi(6,k),nu(6,1),gamma(6,1),phi_r(6,:),phi_p(6,:)] = phi_dp1(q(:,nt),q_d(:,nt),i,j,ax_ib,ay_jb,thetadev);%ar_jb
        end
        %Euler Parameter Constraints
        [phi(7,k),nu(7,1),gamma(7,1),phi_r(7,:),phi_p(7,:)] = phi_ep(q(:,nt),q_d(:,nt),i);
        phi_q = [phi_r phi_p]; %Combined Jacobian
        
        dum1 = 1;
        dum2 = 1;
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
            q(:,nt) = q(:,nt)-(phi_q^-1)*phi(:,k); %Newton-Raphson
            q_d(:,nt) = inv(phi_q)*nu; %Update q_dot
            q_dd(:,nt) = inv(phi_q)*gamma; %update q_doubledot
            k = k+1;
        end
    end
    
%% "ma=F" for initial q_dd values

%Expanded "M" matrix (Only Jp is unknown)
G1 = gmat(q(:,nt),i);
J1p = 4*G1'*J1_b*G1;

%RHS
Bi = bmat(q(:,nt),i,sq_ib);
Bi_d = bmat(q_d(:,nt),i,sq_ib);
A1 = orient(q(:,nt),i);
Fa = [0 0 0]'; %applied forces
F = Fm+Fa; %total external forces
F_b = A1'*F; %In local frame
G1_d = gmat(q_d(:,nt),i);
na_b = [0 0 0]'; %applied torques
nm_b = [0 0 0]'; %distributed torques
n_b = na_b+nm_b; %total external torques
tau_b = 2*G1'*n_b+8*G1_d'*J1_b*G1_d*q([4*(i-1)+3*nb+1:4*(i-1)+3*nb+4],nt); %with no applied torques

%Solve for q_dd and lamda
[q_dd(:,nt),lam(:,nt)] = maf(M,J1p,phi_q,F,tau_b,gamma);
%Reaction forces @ time t=0
R = -phi_q'*lam(:,nt);
Fr_1Op(:,nt) = R(1:3); %Reaction forces at local origin (global reference frame)
nr_Op(:,nt) = A1*0.5*G1*R(4:end); %reaction torques at local origin in local frame
Fr_Q(:,nt) = Fr_1Op(:,nt); %Reaction forces at point Q in global frame
nr_Q(:,nt) = A1*(skew3(-sq_ib)*A1'*Fr_1Op(:,nt))+nr_Op(:,nt); %Reaction torques at point Q in global frame
% clear R Fr_1b nr_1b
%Time Response @ t=0
Op(:,nt) = q(1:3,nt); %O' has coordinates r_i
Op_d(:,nt) = q_d(1:3,nt);
Op_dd(:,nt) = q_dd(1:3,nt);
Q(:,nt) = Op(:,nt)+A1*sq_ib; %Location of end of rod (should be origin)
Q_d(:,nt) = Op_d(:,nt)+Bi*q_d(4:7,nt);
Q_dd(:,nt) = Op_dd(:,nt)+Bi_d*q_d(4:7,nt)+Bi*q_dd(4:7,nt);

%% Dynamic Analysis
for nt = 2:T/h+1
    t(nt) = (nt-1)*h;
    theta(1,nt) = pi/4*cos(2*t(nt))-pi/2;
    theta(2,nt) = -pi/2*sin(2*t(nt));
    theta(3,nt) = -pi*cos(2*t(nt));
    
    %Initial Guess for time nt
    q_dd(:,nt) = q_dd(:,nt-1);
    lam(:,nt) = lam(:,nt-1);
    
    n = 1;
    k = 1;
    clear phi phi_q phi_r phi_p
    while n == 1
        
        %BDF to get q_d and q
        clear r_dd p_dd
        if nt < 3 %not enough for BDF 2
            C = 1;
            beta = 1;
            q_d(:,nt) = C*q_d(:,nt-1)+beta*h*q_dd(:,nt);
            q(:,nt) = C*q(:,nt-1)+beta*h*q_d(:,nt);
        else %BDF 2
            C1 = 4/3; %1 point back
            C2 = -1/3; %2 points back
            beta = 2/3;
            q_d(:,nt) = C1*q_d(:,nt-1)+C2*q_d(:,nt-2)+beta*h*q_dd(:,nt);
            q(:,nt) = C1*q(:,nt-1)+C2*q(:,nt-2)+beta*h*q_d(:,nt);
        end
        
        %% Calculate g(z)

        %Kinematic Constraints
        [phi(1,k),nu(1,1),gamma(1,1),phi_r(1,:),phi_p(1,:)] = phi_dp1(q(:,nt),q_d(:,nt),i,j,az_ib,ay_jb,th2);
        [phi(2,k),nu(2,1),gamma(2,1),phi_r(2,:),phi_p(2,:)] = phi_cd(q(:,nt),q_d(:,nt),i,j,sq_ib,so_jb,c1,dis); %fixed in x
        [phi(3,k),nu(3,1),gamma(3,1),phi_r(3,:),phi_p(3,:)] = phi_cd(q(:,nt),q_d(:,nt),i,j,sq_ib,so_jb,c2,dis); %fixed in y
        [phi(4,k),nu(4,1),gamma(4,1),phi_r(4,:),phi_p(4,:)] = phi_cd(q(:,nt),q_d(:,nt),i,j,sq_ib,so_jb,c3,dis); %fixed in z
        [phi(5,k),nu(5,1),gamma(5,1),phi_r(5,:),phi_p(5,:)] = phi_dp1(q(:,nt),q_d(:,nt),i,j,az_ib,az_jb,th2);
        %Driving Constraints
        if abs(theta(1,nt))*180/pi >= 10 %Change reference angle to avoid singularities
            [phi(6,k),nu(6,1),gamma(6,1),phi_r(6,:),phi_p(6,:)] = phi_dp1(q(:,nt),q_d(:,nt),i,j,ax_ib,ay_jb,theta(:,nt));
        else
            thetadev = theta(:,nt);
            thetadev(1) = theta(1,nt)+pi/4;
            [phi(6,k),nu(6,1),gamma(6,1),phi_r(6,:),phi_p(6,:)] = phi_dp1(q(:,nt),q_d(:,nt),i,j,ax_ib,ay_jb,thetadev);
        end
        %Euler Parameter Constraints
        [phi(7,k),nu(7,1),gamma(7,1),phi_r(7,:),phi_p(7,:)] = phi_ep(q(:,nt),q_d(:,nt),i);
        phi_q = [phi_r phi_p]; %Combined Jacobian
        
        %Evaluate Forces and J Matrix
        G1 = gmat(q(:,nt),i);
        J1p = 4*G1'*J1_b*G1;
        
        Bi = bmat(q(:,nt),i,sq_ib);
        Bi_d = bmat(q_d(:,nt),i,sq_ib);
        A1 = orient(q(:,nt),i);
        Fa = [0 0 0]'; %applied forces
        F = Fm+Fa;
        F_b = A1'*F; %total external forces
        G1_d = gmat(q_d(:,nt),i);
        na_b = [0 0 0]'; %applied torques
        nm_b = [0 0 0]'; %distributed torques
        n_b = na_b+nm_b; %total external torques (about cg)
        tau_b = 2*G1'*n_b+8*G1_d'*J1_b*G1_d*q([4*(i-1)+3*nb+1:4*(i-1)+3*nb+4],nt); %with no applied torques
        
        %Calculate residual (g(z))
        res(:,k) = bdfg(M,J1p,q_dd(:,nt),F,tau_b,h,beta,phi(:,k),phi_q,lam(:,nt)); %calculates residual
        
        %% Convergeance Logic
        
        dum1 = 1;
        dum2 = 1;
        for n2 = 1:size(res,1);
            if abs(res(n2,k)) >= 1e-8 %are any phi values still not zero?
                dum1 = 0; %if so, do not stop while loop
            end
            if k > 1
                if abs(res(n2,k)-res(n2,k-1)) >= 1e-10 %are any phi values still changing significantly?
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
        elseif k > 1000 %if iteration limit is reached
            error('The solver failed to converge')
        else %if not converged
            %Compute psi
            psi = qnpsi(M,J1p,phi_q);
            [r_dd,p_dd] = q2rp(q_dd(:,nt));
            z = [r_dd;p_dd;lam(:,nt)];
            z = z-inv(psi)*res(:,k);
            [q_dd(:,nt),lam(:,nt)] = decz(z,nb);
            clear r_dd p_dd
        end
    end %end of while loop
    %% Desired Quantities @ step nt
    %Reaction Forces/Torques
    R = -phi_q'*lam(:,nt);
    Fr_1Op(:,nt) = R(1:3); %Reaction forces at local origin (global reference frame)
    nr_Op(:,nt) = A1*0.5*G1*R(4:end); %reaction torques at local origin in local frame
    Fr_Q(:,nt) = Fr_1Op(:,nt); %Reaction forces at point Q in global frame
    nr_Q(:,nt) = A1*(skew3(-sq_ib)*A1'*Fr_1Op(:,nt))+nr_Op(:,nt); %Reaction torques at point Q in global frame
    %clear R Fr_1b nr_1b
    
    %Time Response
    Op(:,nt) = q(1:3,nt); %O' has coordinates r_i
    Op_d(:,nt) = q_d(1:3,nt);
    Op_dd(:,nt) = q_dd(1:3,nt);
    Q(:,nt) = Op(:,nt)+A1*sq_ib; %Location of end of rod (should be origin)
    Q_d(:,nt) = Op_d(:,nt)+Bi*q_d(4:7,nt);
    Q_dd(:,nt) = Op_dd(:,nt)+Bi_d*q_d(4:7,nt)+Bi*q_dd(4:7,nt);
end

trun = toc;

%% Plot Results

figure
hold on
plot(t,Op(1,:))
plot(t,Op(2,:))
plot(t,Op(3,:))
title('Position of Point O-prime')
ylabel('Position [m]')
xlabel('Time [s]')
legend('Global X','Global Y','Global Z')

figure
hold on
plot(t,Op_d(1,:))
plot(t,Op_d(2,:))
plot(t,Op_d(3,:))
title('Velocity of Point O-prime')
ylabel('Velocity [m/s]')
xlabel('Time [s]')
legend('Global X','Global Y','Global Z')

figure
hold on
plot(t,Op_dd(1,:))
plot(t,Op_dd(2,:))
plot(t,Op_dd(3,:))
title('Acceleration of Point O-prime')
ylabel('Acceleration [m/s^2]')
xlabel('Time [s]')
legend('Global X','Global Y','Global Z')


figure
hold on
plot(t,Q(1,:))
plot(t,Q(2,:))
plot(t,Q(3,:))
title('Position of Point Q')
ylabel('Position [m]')
xlabel('Time [s]')
legend('Global X','Global Y','Global Z')

figure
hold on
plot(t,Q_d(1,:))
plot(t,Q_d(2,:))
plot(t,Q_d(3,:))
title('Velocity of Point Q')
ylabel('Velocity [m/s]')
xlabel('Time [s]')
legend('Global X','Global Y','Global Z')

figure
hold on
plot(t,Q_dd(1,:))
plot(t,Q_dd(2,:))
plot(t,Q_dd(3,:))
title('Acceleration of Point Q')
ylabel('Acceleration [m/s^2]')
xlabel('Time [s]')
legend('Global X','Global Y','Global Z')

figure
hold on
plot(t,Fr_Q(1,:))
plot(t,Fr_Q(2,:))
plot(t,Fr_Q(3,:))
title('Reaction Forces')
ylabel('Forces [N]')
xlabel('Time [s]')
legend('Global X','Global Y','Global Z')

figure
hold on
plot(t,nr_Q(1,:))
plot(t,nr_Q(2,:))
plot(t,nr_Q(3,:))
title(['Reaction Torques (h = ' num2str(h) ')'])
ylabel('Torques [N-m]')
xlabel('Time [s]')
legend('Global X','Global Y','Global Z')