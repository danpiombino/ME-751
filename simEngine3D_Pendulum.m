% simEngine3D_Pendulum.m
% creates model for HW6 Problem 3
%
% Written by: Dan Piombino
%12/6/16

clear;
clc;
tic;

%% Constants, Etc.
%Environment
g = 9.81*[0 -1 0]'; %gravity in m/s^2

%System
method = 3;
nb = 1; %number of links
nci = 6; %number of constraints including initial conditions
nc = 5; %number of constraints after initial conditions
T = 2e3; %2e3 ref
nstep = 4e3*20; %1e5; %4e3 ref (must be multiple to compare solutions)
h = T/nstep; %time step

%Reference Solution
Ref = csvread('Simple Pendulum Solution.csv');
tref = Ref(:,1);
xref = Ref(:,2);
yref = Ref(:,3);
href = tref(2)-tref(1);
clear Ref

%Define bodies (Mass and geometric properties)
L = 1; %1m rod length
m = 1; %1kg lumped mass
M_body = m*eye(3); %rod mass matrix
J_b = zeros(3,3);

%"Guess" at Initial Positions/Orientations
q(:,1) = -L*[1 0 0]';
q([4:7],1) = [0 0 0 1]';
q_d = zeros(7*nb,1);
q_dd = zeros(7*nb,1);

%Define Combined Mass Properties (Time Invariant)
M = zeros(3*nb,3*nb);
for i = 1:nb
    M(3*(i-1)+1:3*(i-1)+3,3*(i-1)+1:3*(i-1)+3) = M_body;
end

%Define Time Invariant Mass-Distributed + Applied Forces
F_c = zeros(3*nb,1); %for speed
for i = 1:nb
    F_c(3*(i-1)+1:3*(i-1)+3,1) = m*g; %constant forces
end

%Define Time Invariant Applied Torques (omega form)
n_bc = zeros(3*nb,1); %no applied torques

%Reference Vectors
a0 = [0 0 0]';
ax = [1 0 0]';
ay = [0 1 0]';
az = [0 0 1]';
dis_0 = [0 0 0]'; %reference distance for cd's
dis_i = [0 -1 0]'; %initial distance, velocity, acceleration of B0 in x direction
th_c = [1 0 0]'; %unscaled constant theta for constraint input (multiply by desired angle in rad's)

%% Initial Constraints

nt = 1; %time step number
t(nt) = 0;
dum = 1;
iter(nt) = 0; %iteration number
while dum == 1
    iter(nt) = iter(nt)+1;
    n = 1; %constraint number
    %Initial Condition Constraint
    [phi(n,1),nu(n,1),gamma(n,1),phi_r(n,:),phi_p(n,:),~] = phi_cd(q(:,nt),q_d(:,nt),1,0,dis_0,-L*ax,ax,dis_0,[]); %x location to ground
    n = n+1;
    %CD's to ground
    [phi(n,1),nu(n,1),gamma(n,1),phi_r(n,:),phi_p(n,:),~] = phi_cd(q(:,nt),q_d(:,nt),1,0,-L*ax,dis_0,ax,dis_0,[]); %x location to ground
    n = n+1;
    [phi(n,1),nu(n,1),gamma(n,1),phi_r(n,:),phi_p(n,:),~] = phi_cd(q(:,nt),q_d(:,nt),1,0,-L*ax,dis_0,ay,dis_0,[]); %y location to ground
    n = n+1;
    [phi(n,1),nu(n,1),gamma(n,1),phi_r(n,:),phi_p(n,:),~] = phi_cd(q(:,nt),q_d(:,nt),1,0,-L*ax,dis_0,az,dis_0,[]); %z location to ground
    n = n+1;
    %DP1's to ground
    [phi(n,1),nu(n,1),gamma(n,1),phi_r(n,:),phi_p(n,:),~] = phi_dp1(q(:,nt),q_d(:,nt),1,0,az,ax,pi/2*th_c,[]); %local z perpendicular to global x
    n = n+1;
    [phi(n,1),nu(n,1),gamma(n,1),phi_r(n,:),phi_p(n,:),~] = phi_dp1(q(:,nt),q_d(:,nt),1,0,az,ay,pi/2*th_c,[]); %local z perpendicular to global y
    n = n+1;
    for i = 1:nb
        [phi(n,1),nu(n,1),gamma(n,1),phi_r(n,:),phi_p(n,:),~] = phi_ep(q(:,nt),q_d(:,nt),i,[]); %Euler parameterization constraints
        n = n+1;
    end
    n = n-1;
    phi_q = [phi_r phi_p];
    
    %Convergence Logic
    if max(abs(phi)) <= 1e-5 %consistent position IC's
        if max(abs(phi_q*q_d(:,nt)-nu)) <= 1e-5 %consistent velocity IC's
            dum = 0; %Fully consistent IC's
        else
            q_d(:,nt) = inv(phi_q)*nu;
        end
    else
        q(:,nt) = q(:,nt)-inv(phi_q)*phi; %Newton-Raphson
        q_d(:,nt) = inv(phi_q)*nu; %Update q_dot
    end
end

%% Initial Configuration Inverse Dynamics
%Remove IC Constraint
phi = phi(2:end);
phi_q = phi_q(2:end,:);
phi_r = phi_r(2:end,:);
phi_p = phi_p(2:end,:);
nu = nu(2:end);
gamma = gamma(2:end);

%Applied Forces/Torques (Time Dependent)
F_t = zeros(3*nb,1); %no time dependent forces
n_bt = zeros(3*nb,1); %no time dependent torques
F = F_c + F_t; %total applied forces at current time step
n_b = n_bc + n_bt; %total applied torques at current time step

%Inertia Matrix
Jp = zeros(4*nb,4*nb);
tau_b = zeros(4*nb,1);
for i = 1:nb
    G = gmat(q(:,nt),i);
    G_d = gmat(q_d(:,nt),i);
    Jp(4*(i-1)+1:4*(i-1)+4,4*(i-1)+1:4*(i-1)+4) = 4*G'*J_b*G;
    tau_b(4*(i-1)+1:4*(i-1)+4,1) = 2*G'*n_b(3*(i-1)+1:3*(i-1)+3,1)+8*G_d'*J_b*G_d*q(3*nb+4*(i-1)+1:3*nb+4*(i-1)+4,nt);
    clear G G_d
end

%Solve inverse dynamics for accelerations and lambdas
[q_dd(:,nt),lam(:,nt)] = maf(M,Jp,phi_q,F,tau_b,gamma(:,nt));
F_r(:,nt) = -phi_r'*lam(:,nt); %Reaction forces @ current time
tau_br(:,nt) = -phi_p'*lam(:,nt); %Reaction torques @ current time (r-p)

%Desired quantities at t=0
mt = 1; %reference solution vector index
dx(mt) = q(1,nt)-xref(mt);
dy(mt) = q(2,nt)-yref(mt);
mt = mt+1;

%% Dynamic Analysis
for nt = 2:nstep+1
    clc
    nt
    t(nt) = (nt-1)*h;
    
    %Initial Guess for time nt
    q_dd(:,nt) = q_dd(:,nt-1);
    lam(:,nt) = lam(:,nt-1);
    
    dum = 1;
    k = 1;
    clear res norm
    while dum == 1
        
        %BDF to get q_d and q
        clear r_dd p_dd phi phi_q phi_r phi_p nu gamma
        if nt > 7 %BDF 6
            C1 = 360/147;
            C2 = -450/147;
            C3 = 400/147;
            C4 = -225/147;
            C5 = 72/147;
            C6 = -10/147;
            beta = 60/147;
            q_d(:,nt) = C1*q_d(:,nt-1)+C2*q_d(:,nt-2)+C3*q_d(:,nt-3)+C4*q_d(:,nt-4)+C5*q_d(:,nt-5)+C6*q_d(:,nt-6)+beta*h*q_dd(:,nt);
            q(:,nt) = C1*q(:,nt-1)+C2*q(:,nt-2)+C3*q(:,nt-3)+C4*q(:,nt-4)+C5*q(:,nt-5)+C6*q(:,nt-6)+beta*h*q_d(:,nt);
        elseif nt > 6 %BDF 5
%         if nt > 6 %BDF 5
            C1 = 300/137;
            C2 = -300/137;
            C3 = 200/137;
            C4 = -75/137;
            C5 = 12/137;
            beta = 60/137;
            q_d(:,nt) = C1*q_d(:,nt-1)+C2*q_d(:,nt-2)+C3*q_d(:,nt-3)+C4*q_d(:,nt-4)+C5*q_d(:,nt-5)+beta*h*q_dd(:,nt);
            q(:,nt) = C1*q(:,nt-1)+C2*q(:,nt-2)+C3*q(:,nt-3)+C4*q(:,nt-4)+C5*q(:,nt-5)+beta*h*q_d(:,nt);
        elseif nt > 5 %BDF 4
%         if nt > 5 %BDF 4
            C1 = 48/25;
            C2 = -36/25;
            C3 = 16/25;
            C4 = -3/25;
            beta = 12/25;
            q_d(:,nt) = C1*q_d(:,nt-1)+C2*q_d(:,nt-2)+C3*q_d(:,nt-3)+C4*q_d(:,nt-4)+beta*h*q_dd(:,nt);
            q(:,nt) = C1*q(:,nt-1)+C2*q(:,nt-2)+C3*q(:,nt-3)+C4*q(:,nt-4)+beta*h*q_d(:,nt);
        elseif nt > 4 %BDF 3
%         if nt > 4 %BDF 3
            C1 = 18/11;
            C2 = -9/11;
            C3 = 2/11;
            beta = 6/11;
            q_d(:,nt) = C1*q_d(:,nt-1)+C2*q_d(:,nt-2)+C3*q_d(:,nt-3)+beta*h*q_dd(:,nt);
            q(:,nt) = C1*q(:,nt-1)+C2*q(:,nt-2)+C3*q(:,nt-3)+beta*h*q_d(:,nt);
        elseif nt > 3 %BDF 2
%         if nt > 3
            C1 = 4/3; %1 point back
            C2 = -1/3; %2 points back
            beta = 2/3;
            q_d(:,nt) = C1*q_d(:,nt-1)+C2*q_d(:,nt-2)+beta*h*q_dd(:,nt);
            q(:,nt) = C1*q(:,nt-1)+C2*q(:,nt-2)+beta*h*q_d(:,nt);
        else %BDF 1 (Backward Euler)
            C = 1;
            beta = 1;
            q_d(:,nt) = C*q_d(:,nt-1)+beta*h*q_dd(:,nt);
            q(:,nt) = C*q(:,nt-1)+beta*h*q_d(:,nt);
        end
        
        %% Calculate g(z)
        
        clear phi phi_q phi_r phi_p gamma nu rf_q rf_p
        %Kinematic Constraints
        n = 1;
        %CD's to ground
        [phi(n,1),nu(n,1),gamma(n,1),phi_r(n,:),phi_p(n,:),rf_q(:,:,n)] = phi_cd(q(:,nt),q_d(:,nt),1,0,-L*ax,dis_0,ax,dis_0,lam(n,nt)); %x location to ground
        n = n+1;
        [phi(n,1),nu(n,1),gamma(n,1),phi_r(n,:),phi_p(n,:),rf_q(:,:,n)] = phi_cd(q(:,nt),q_d(:,nt),1,0,-L*ax,dis_0,ay,dis_0,lam(n,nt)); %y location to ground
        n = n+1;
        [phi(n,1),nu(n,1),gamma(n,1),phi_r(n,:),phi_p(n,:),rf_q(:,:,n)] = phi_cd(q(:,nt),q_d(:,nt),1,0,-L*ax,dis_0,az,dis_0,lam(n,nt)); %z location to ground
        n = n+1;
        %DP1's to ground
        [phi(n,1),nu(n,1),gamma(n,1),phi_r(n,:),phi_p(n,:),rf_q(:,:,n)] = phi_dp1(q(:,nt),q_d(:,nt),1,0,az,ax,pi/2*th_c,lam(n,nt)); %local z perpendicular to global x
        n = n+1;
        [phi(n,1),nu(n,1),gamma(n,1),phi_r(n,:),phi_p(n,:),rf_q(:,:,n)] = phi_dp1(q(:,nt),q_d(:,nt),1,0,az,ay,pi/2*th_c,lam(n,nt)); %local z perpendicular to global y
        n = n+1;
        %EP Normalization Constraints
        for i = 1:nb
            [phi(n,1),nu(n,1),gamma(n,1),phi_r(n,:),phi_p(n,:),rf_p(:,:,i)] = phi_ep(q(:,nt),q_d(:,nt),i,lam(n,nt)); %Euler parameterization constraints
            n = n+1;
        end
        n = n-1;
        %Combined Constraint Jacobian
        phi_q = [phi_r phi_p];
        
        %Applied Forces/Torques (Time Dependent)
        F_t = zeros(3*nb,1); %no time dependent forces
        n_bt = zeros(3*nb,1); %no time dependent torques
        F = F_c + F_t; %total applied forces at current time step
        n_b = n_bc + n_bt; %total applied torques at current time step
        
        %Inertia Matrix
        Jp = zeros(4*nb,4*nb);
        tau_b = zeros(4*nb,1);
        for i = 1:nb
            G = gmat(q(:,nt),i);
            G_d = gmat(q_d(:,nt),i);
            Jp(4*(i-1)+1:4*(i-1)+4,4*(i-1)+1:4*(i-1)+4) = 4*G'*J_b*G;
            tau_b(4*(i-1)+1:4*(i-1)+4,1) = 2*G'*n_b(3*(i-1)+1:3*(i-1)+3,1)+8*G_d'*J_b*G_d*q(3*nb+4*(i-1)+1:3*nb+4*(i-1)+4,nt);
            clear G G_d
        end
        
        %Method Dependet Calculations
        if method == 1
%             Full NR info (Valid for no applied forces)
            %Partials of Force
            F_r = zeros(3*nb,3*nb); %for constant forces
            F_rd = zeros(3*nb,3*nb); %for constant forces
            F_p = zeros(3*nb,4*nb); %for constant forces
            F_pd = zeros(3*nb,4*nb); %for constant forces
            %Partials of Torques
            tau_b_r = zeros(4*nb,3*nb); %for constant torques
            tau_b_rd = zeros(4*nb,3*nb); %for constant torques
            
            tau_b_p = zeros(4*nb,4*nb); %for speed
            tau_b_pd = zeros(4*nb,4*nb); %for speed
            Jp_p = zeros(4*nb,4*nb);
            for i = 1:nb
                G = gmat(q(:,nt),i);
                G_d = gmat(q_d(:,nt),i);
                G_dd = gmat(q_dd(:,nt),i);
                alpha = J_b*G_d*q([4*(1-1)+3*nb+1:4*(1-1)+3*nb+4],nt);
                tau_b_p(4*(i-1)+1:4*(i-1)+4,4*(i-1)+1:4*(i-1)+4) = 8*G_d'*J_b*G_d; %for only inertial torques
                tau_b_pd(4*(i-1)+1:4*(i-1)+4,4*(i-1)+1:4*(i-1)+4) = 8*tmat(alpha); %for only inertial torques
                Jp_p(4*(i-1)+1:4*(i-1)+4,4*(i-1)+1:4*(i-1)+4) = -4*G'*J_b*G_dd+4*tmat(alpha); %for only inertial torques
            end
            
            %Combining partials of reaction forces
            rf_q_full = zeros(size(rf_q,1),size(rf_q,2)); %reaction forces associated with basic constraints
            for jn = 1:nc
                rf_q_full = rf_q_full + rf_q(:,:,jn);
            end
            rf_p_full = zeros(size(rf_p,1),size(rf_p,2)); %reaction forces associated with EP constraints
            for jn = 1:nb
                rf_p_full = rf_p_full+rf_p(:,:,jn);
            end
            
            %Calculate psi
            psi = nrpsi(M,h,beta,Jp,Jp_p,phi_q,F_r,F_rd,F_p,F_pd,tau_b_r,tau_b_rd,tau_b_p,tau_b_pd,rf_q_full,rf_p_full);
        elseif method == 2
            if k == 1 %must keep track of iteration # for each time step
                %Full NR info (Valid for no applied forces)
                %Partials of Force
                F_r = zeros(3*nb,3*nb); %for constant forces
                F_rd = zeros(3*nb,3*nb); %for constant forces
                F_p = zeros(3*nb,4*nb); %for constant forces
                F_pd = zeros(3*nb,4*nb); %for constant forces
                %Partials of Torques
                tau_b_r = zeros(4*nb,3*nb); %for constant torques
                tau_b_rd = zeros(4*nb,3*nb); %for constant torques
                
                tau_b_p = zeros(4*nb,4*nb); %for speed
                tau_b_pd = zeros(4*nb,4*nb); %for speed
                Jp_p = zeros(4*nb,4*nb);
                for i = 1:nb
                    G = gmat(q(:,nt),i);
                    G_d = gmat(q_d(:,nt),i);
                    G_dd = gmat(q_dd(:,nt),i);
                    alpha = J_b*G_d*q([4*(1-1)+3*nb+1:4*(1-1)+3*nb+4],nt);
                    tau_b_p(4*(i-1)+1:4*(i-1)+4,4*(i-1)+1:4*(i-1)+4) = 8*G_d'*J_b*G_d; %for only inertial torques
                    tau_b_pd(4*(i-1)+1:4*(i-1)+4,4*(i-1)+1:4*(i-1)+4) = 8*tmat(alpha); %for only inertial torques
                    Jp_p(4*(i-1)+1:4*(i-1)+4,4*(i-1)+1:4*(i-1)+4) = -4*G'*J_b*G_dd+4*tmat(alpha); %for only inertial torques
                end
                
                %Combining partials of reaction forces
                rf_q_full = zeros(size(rf_q,1),size(rf_q,2)); %reaction forces associated with basic constraints
                for jn = 1:nc
                    rf_q_full = rf_q_full + rf_q(:,:,jn);
                end
                rf_p_full = zeros(size(rf_p,1),size(rf_p,2)); %reaction forces associated with EP constraints
                for jn = 1:nb
                    rf_p_full = rf_p_full+rf_p(:,:,jn);
                end
                
                %Calculate psi
                psi = nrpsi(M,h,beta,Jp,Jp_p,phi_q,F_r,F_rd,F_p,F_pd,tau_b_r,tau_b_rd,tau_b_p,tau_b_pd,rf_q_full,rf_p_full);
            end
        else
            %Quasi-Newton
            psi = qnpsi(M,Jp,phi_q);
        end
        
        
        %Calculate residual (g(z))
        res(:,1) = bdfg(M,Jp,q_dd(:,nt),F,tau_b,h,beta,phi(:,1),phi_q,lam(:,nt)); %calculates residual
        dz = -inv(psi)*res(:,1);
        
        %% Convergeance Logic
        dum1 = 1;
        if nt > 6 %Higher order BDF methods
            if abs(norm(dz)) >= 1e-3 || abs(norm(phi)) >= 1e-3 || abs(norm(phi_q*q_d(:,nt)-nu)) >= 1e-3 %are any phi values still not zero?
                dum1 = 0; %if so, do not stop while loop
            end
        else
            if abs(norm(dz)) >= 1e-3 || abs(norm(phi)) >= 1e-3 %are any phi values still not zero?
                dum1 = 0; %if so, do not stop while loop
            end
        end
        if dum1 == 1 %if every value of res is below threshold
            dum = 0;
            iter(nt) = k;
        elseif k > 1000 %if iteration limit is reached
            error('The solver failed to converge')
        else %if not converged
            z = [q_dd(:,nt);lam(:,nt)];
            z = z+dz;
            [q_dd(:,nt),lam(:,nt)] = decz(z,nb);
            k = k+1;
            clear r_dd p_dd
        end
    end %end of while loop
    %% Desired Quantities @ step nt
    
    %Velocity Constraints
    vc(:,nt) = phi_q*q_d(:,nt)-nu;
    nvc(nt) = norm(vc(:,nt));
    
    %Reaction forces & Torques
    F_r(:,nt) = -phi_r'*lam(:,nt); %Reaction forces @ current time
    tau_br(:,nt) = -phi_p'*lam(:,nt); %Reaction torques @ current time (r-p)
    
    %Check with Reference Solution
    if abs(round(t(nt)/href) - t(nt)/href) < h/10
        dx(mt) = q(1,nt)-xref(mt);
        dy(mt) = q(2,nt)-yref(mt);
        mt = mt+1;
    end
    
end

% if max(abs(dx)) || max(abs(dy)) >= 1e-1 %low precision tolerance
%     disp('Low Precision: Failure')
%     disp('High Precision: Failure')
% elseif max(abs(dx)) || max(abs(dy)) >= 1e-3
%     disp('Low Precision: Pass')
%     disp('High Precision: Failure')
% else
%     disp('Low Precision: Pass')
%     disp('High Precision: Pass')
% end

mdx = max(abs(dx));
adx = mean(abs(dx));
mdy = max(abs(dy));
ady = mean(abs(dy));
mit = max(iter);
ait = mean(iter);
Tsi = 2.365;
Tsf = 2.3607;
Tri = 2.3596;
Trf = 2.3546;

trun = toc
cd('C:\Users\Dan\Documents\Classwork\ME-751\SimEngine3D\Final Project\Final Results');
% save('Pendulum_QN');

%% Plot Results

figure
hold on
plot(t([1:10/T*nstep+1]),q(1,[1:10/T*nstep+1]),'k')
plot(tref(1:10/2000*4e3+1),xref(1:10/2000*4e3+1),'k.')
plot(t([1:10/T*nstep+1]),q(2,[1:10/T*nstep+1]),'b')
plot(tref(1:10/2000*4e3+1),yref(1:10/2000*4e3+1),'b.')
% title('First 10s Comparison')
ylabel('Lumped Mass Position [m]')
xlabel('Time [s]')
legend('Sim X','Ref X','Sim Y','Ref Y')

figure
hold on
plot(t([end-10/T*nstep:end]),q(1,[end-10/T*nstep:end]),'k','LineWidth',1)
plot(tref(end-10/2000*4e3:end),xref(end-10/2000*4e3:end),'k-.')
plot(t([end-10/T*nstep:end]),q(2,[end-10/T*nstep:end]),'b','LineWidth',1)
plot(tref(end-10/2000*4e3:end),yref(end-10/2000*4e3:end),'b-.')
xlim([T-10 T]);
% title('Last 10s Comparison')
ylabel('Lumped Mass Position [m]')
xlabel('Time [s]')
legend('Sim X','Ref X','Sim Y','Ref Y')

figure
hold on
plot(t([1:10/T*nstep+1]),q(1,[1:10/T*nstep+1]),'k')
plot(tref(1:10/2000*4e3+1),xref(1:10/2000*4e3+1),'k.')
plot(t([1:10/T*nstep+1]),q(1,[end-10/T*nstep:end]),'b')
plot(tref(1:10/2000*4e3+1),xref(end-10/2000*4e3:end),'b-.')
% title('Initial and Final 10s Comparison')
ylabel('Lumped Mass X Position [m]')
xlabel('Time Into Window [s]')
legend('Initial Sim','Initial Ref','Final Sim','Final Ref')

figure %may or may not be meaningful
hold on
plot(tref,dx)
plot(tref,dy)
% title('Simulation Error Compared to Reference Solution')
ylabel('Position Difference [m]')
xlabel('Time [s]')
legend('X Error','Y Error')

% figure %simulated positions (start)
% hold on
% plot(t([1:10/T*nstep+1]),q(1,[1:10/T*nstep+1]))
% plot(t([1:10/T*nstep+1]),q(2,[1:10/T*nstep+1]))
% plot(t([1:10/T*nstep+1]),q(3,[1:10/T*nstep+1]))
% title('Position of Point P1 (First 10s)')
% ylabel('Position [m]')
% xlabel('Time [s]')
% legend('Global X','Global Y','Global Z')
% 
% figure %simulated positions (end)
% hold on
% plot(t([end-10/T*nstep:end]),q(1,[end-10/T*nstep:end]))
% plot(t([end-10/T*nstep:end]),q(2,[end-10/T*nstep:end]))
% plot(t([end-10/T*nstep:end]),q(3,[end-10/T*nstep:end]))
% title('Position of Point P1 (Last 10s)')
% ylabel('Position [m]')
% xlabel('Time [s]')
% legend('Global X','Global Y','Global Z')
% xlim([1990 2000]);



% figure
% hold on
% plot(t([end-10/T*nstep:end]),q(1,[end-10/T*nstep:end]))
% plot(tref(end-10/T*4e3:end),xref(end-10/T*4e3:end),'--')
% title('Global X of Point P1 (Last 10s)')
% ylabel('Position [m]')
% xlabel('Time [s]')
% legend('simEngine3D','Reference')
% xlim([1990 2000]);
% 
% figure
% hold on
% plot(t([end-10/T*nstep:end]),q(2,[end-10/T*nstep:end]))
% plot(tref(end-10/T*4e3:end),yref(end-10/T*4e3:end),'--')
% title('Global Y of Point P1 (Last 10s)')
% ylabel('Position [m]')
% xlabel('Time [s]')
% legend('simEngine3D','Reference')
% xlim([1990 2000]);

%Compare New vs old
% figure
% hold on
% plot(t([1:10/T*nstep+1]),q(1,[1:10/T*nstep+1]),'k')
% plot(t([1:10/T*nstep+1]),q(1,[end-10/T*nstep:end]),'k--')
% plot(tref(1:10/T*4e3+1),xref(1:10/T*4e3+1),'b')
% plot(tref(1:10/T*4e3+1),xref(end-10/T*4e3:end),'b--')
% title('X Position Time Window Comparison')
% ylabel('Position [m]')
% xlabel('Time [s]')
% legend('simEngine3D Start','simEngine3D End','Reference Start','Reference End')



% figure
% hold on
% plot(t,q_d(1,:))
% plot(t,q_d(2,:))
% plot(t,q_d(3,:))
% title('Velocity of Point P1')
% ylabel('Velocity [m/s]')
% xlabel('Time [s]')
% legend('Global X','Global Y','Global Z')

% figure
% hold on
% plot(t,nvc)
% title('Velocity Constraint Violation')
% ylabel('Norm of Violations')
% xlabel('Time [s]')
% 
% figure
% plot(t,iter)
% % title('Iterations to Convergence')
% ylabel('Number of Iterations')
% xlabel('Time [s]')