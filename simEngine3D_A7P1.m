% simEngine3D_A7P1.m
% creates model for HW6 Problem 3
%
% Written by: Dan Piombino
% 10/27/16

clear;
clc;

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
dt = 1e-3; %time step
T = 10; %length of simulation
nb = 1;

%% Mass and other properties

rho = 7800;
w = .05; %x-section width
A = w^2; %area
m = A*L*rho;
M = eye(3)*m; %since there's only 1 body
Fa = m*9.81*[0 0 -1]';

J1_bx = m/12*(w^2+w^2);
J1_by = m/12*(w^2+(2*L)^2);
J1_bz = J1_by;
J1_b = [J1_bx 0 0; 0 J1_by 0; 0 0 J1_bz];

%% Kinematic Analysis
for nt = 1:T/dt+1
    t(nt) = (nt-1)*dt;
    theta(1,nt) = pi/4*cos(2*t(nt));
    theta(2,nt) = -pi/2*sin(2*t(nt));
    theta(3,nt) = -pi*cos(2*t(nt));
    
    n = 1;
    k = 1;
    clear phi phi_q phi_r phi_p
    while n == 1
        %Kinematic Constraints
        [phi(1,k),nu(1,1),gamma(1,1),phi_r(1,:),phi_p(1,:)] = phi_dp1(q(:,nt),q_d(:,nt),i,j,az_ib,ay_jb,th2);
        [phi(2,k),nu(2,1),gamma(2,1),phi_r(2,:),phi_p(2,:)] = phi_cd(q(:,nt),q_d(:,nt),i,j,sq_ib,so_jb,c1,dis); %fixed in x
        [phi(3,k),nu(3,1),gamma(3,1),phi_r(3,:),phi_p(3,:)] = phi_cd(q(:,nt),q_d(:,nt),i,j,sq_ib,so_jb,c2,dis); %fixed in y
        [phi(4,k),nu(4,1),gamma(4,1),phi_r(4,:),phi_p(4,:)] = phi_cd(q(:,nt),q_d(:,nt),i,j,sq_ib,so_jb,c3,dis); %fixed in z
        [phi(5,k),nu(5,1),gamma(5,1),phi_r(5,:),phi_p(5,:)] = phi_dp1(q(:,nt),q_d(:,nt),i,j,az_ib,az_jb,th2);
        %Driving Constraints
        if abs(theta(1,nt))*180/pi >= 10 %Change reference angle to avoid singularities
            [phi(6,k),nu(6,1),gamma(6,1),phi_r(6,:),phi_p(6,:)] = phi_dp1(q(:,nt),q_d(:,nt),i,j,ax_ib,-az_jb,theta(:,nt));
        else
            thetadev = theta(:,nt);
            thetadev(1) = theta(1,nt)+pi/4;
            [phi(6,k),nu(6,1),gamma(6,1),phi_r(6,:),phi_p(6,:)] = phi_dp1(q(:,nt),q_d(:,nt),i,j,ax_ib,ar_jb,thetadev);
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
            iter(nt) = k;
        elseif k > 1000 %if iteration limit is reached
            error('The solver failed to converge')
        else
            q(:,nt) = q(:,nt)-(phi_q^-1)*phi(:,k); %Newton-Raphson
            q_d(:,nt) = inv(phi_q)*nu; %Update q_dot
            q_dd(:,nt) = inv(phi_q)*gamma; %update q_doubledot
            k = k+1;
        end
    end
    
    %Desired Quantities @ step nt
    Bi = bmat(q(:,nt),i,sq_ib);
    Bi_d = bmat(q_d(:,nt),i,sq_ib);
    A = orient(q(:,nt),i);
    
    G1 = gmat(q(:,nt),i);
    G1_d = gmat(q_d(:,nt),i);
    J1p = 4*G1'*J1_b*G1;
    Fa_b = A'*Fa; %Applied Forces to Body i in local reference frame
    tau_b = 8*G1_d'*J1_b*G1_d*q([4*(i-1)+3*nb+1:4*(i-1)+3*nb+4],nt); %with no applied torques
    
    idynLHS = [phi_r([1:6],:)' zeros(3*nb,nb);phi_p([1:6],:)' phi_p(7,:)']; %inverse dynamics left hand side
    idynRHS = [M*q([3*(i-1)+1:3*(i-1)+3],nt)-Fa_b; J1p*q_dd([4*(i-1)+3*nb+1:4*(i-1)+3*nb+4],nt)-tau_b]; %inverse dynamics right hand side
    lam(:,nt) = -inv(idynLHS)*idynRHS;
    Fr_b(:,nt) = -phi_r'*lam(:,nt);
    nr_Ob(:,nt) = -0.5*G1*phi_p'*lam(:,nt); %converts to r-w
    Fr = A*Fr_b;
    nr_O = A*nr_Ob;
    nr = nr_O+skew3(q(1:3,nt))*Fr; %reaction torque at point Q
    
    Op(:,nt) = q(1:3,nt); %O' has coordinates r_i
    Op_d(:,nt) = q_d(1:3,nt);
    Op_dd(:,nt) = q_dd(1:3,nt);
    Q(:,nt) = Op(:,nt)+A*sq_ib; %Location of end of rod (should be origin)
    Q_d(:,nt) = Op_d(:,nt)+Bi*q_d(4:7,nt);
    Q_dd(:,nt) = Op_dd(:,nt)+Bi_d*q_d(4:7,nt)+Bi*q_dd(4:7,nt);
    
    %Set up next iteration
    if nt < T/dt+1
        q(:,nt+1) = q(:,nt); %initial guess for nt+1 = converged answer for nt
        q_d(:,nt+1) = q_d(:,nt);
        q_dd(:,nt+1) = q_dd(:,nt);
    end
end

%% Plot

% figure
% hold on
% plot(t,Op(1,:))
% plot(t,Op(2,:))
% plot(t,Op(3,:))
% title('Position of Point O-prime')
% ylabel('Position [m]')
% xlabel('Time [s]')
% legend('Global X','Global Y','Global Z')
% 
% figure
% hold on
% plot(t(1:length(Op_d)),Op_d(1,:))
% plot(t(1:length(Op_d)),Op_d(2,:))
% plot(t(1:length(Op_d)),Op_d(3,:))
% title('Velocity of Point O-prime')
% ylabel('Velocity [m/s]')
% xlabel('Time [s]')
% legend('Global X','Global Y','Global Z')
% 
% figure
% hold on
% plot(t(1:length(Op_dd)),Op_dd(1,:))
% plot(t(1:length(Op_dd)),Op_dd(2,:))
% plot(t(1:length(Op_dd)),Op_dd(3,:))
% title('Acceleration of Point O-prime')
% ylabel('Acceleration [m/s^2]')
% xlabel('Time [s]')
% legend('Global X','Global Y','Global Z')
% 
% 
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
% plot(t(1:length(Q_d)),Q_d(1,:))
% plot(t(1:length(Q_d)),Q_d(2,:))
% plot(t(1:length(Q_d)),Q_d(3,:))
% title('Velocity of Point Q')
% ylabel('Velocity [m/s]')
% xlabel('Time [s]')
% legend('Global X','Global Y','Global Z')
% 
% figure
% hold on
% plot(t(1:length(Q_dd)),Q_dd(1,:))
% plot(t(1:length(Q_dd)),Q_dd(2,:))
% plot(t(1:length(Q_dd)),Q_dd(3,:))
% title('Acceleration of Point Q')
% ylabel('Acceleration [m/s^2]')
% xlabel('Time [s]')
% legend('Global X','Global Y','Global Z')

figure
hold on
plot(t,Fr_b(1,:))
plot(t,Fr_b(2,:))
plot(t,Fr_b(3,:))
title('Reaction Forces')
ylabel('Forces [N]')
xlabel('Time [s]')
legend('Local X','Local Y','Local Z')

figure
hold on
plot(t,nr(1,:))
plot(t,nr(2,:))
plot(t,nr(3,:))
title('Reaction Torques')
ylabel('Torques [N-m]')
xlabel('Time [s]')
legend('Local X','Local Y','Local Z')