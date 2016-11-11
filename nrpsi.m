function psi = nrpsi(M,h,beta,Jp,Jp_p,phi_q,F_r,F_rd,F_p,F_pd,tau_r,tau_rd,tau_p,tau_pd,rf_q,rf_p)
% Evaluates full jacobian (psi) for newton-raphson or modified newton
% 
% Inputs:   M: Mass matrix [3*nb x 3*nb]
%           h: Step size
%           beta: BDF method beta constant
%           Jp: Euler Parameter Inertia Matrix [4*nb x 4*nb]
%           phi_q: Jacobian of constraint functions [nc+nb x 7*nb]
%           F_: Partial derivative of Forces with respect to (r,r_d,p,p_d)
%           tau_: Partial derivative of torques with respect to
%                   (r,r_d,p,p_d)
%           rf_q: Partial derivative of reaction forces with respect to q
%           rf_p: Partial derivative of euler parameter reaction forces
%                   with respect to q
%
% Outputs:   psi: Jacobian of g(z) for newton-raphson or modified newton
%                   method
% 
% Written By: Dan Piombino
% 11/10/16

nb = size(M,1)/3;
nc = size(phi_q,1)-nb;

rf_rr = rf_q(1:3*nb,1:3*nb);
rf_rp = rf_q(1:3*nb,3*nb+1:end);
rf_pp = rf_q(3*nb+1:end,3*nb+1:end);
rf_pr = rf_q(3*nb+1:end,1:3*nb);

psi_11 = M+h^2*beta^2*rf_rr-h^2*beta^2*F_r-h*beta*F_rd;
psi_12 = h^2*beta^2*rf_rp-h^2*beta^2*F_p-h*beta*F_pd;
psi_21 = h^2*beta^2*rf_pr-h^2*beta^2*tau_r-h*beta*tau_rd;
psi_22 = Jp+h^2*beta^2*Jp_p+h^2*beta^2*rf_p+h^2*beta^2*rf_pp-h^2*beta^2*tau_p-h*beta*tau_pd;

psi_1 = [psi_11 psi_12;psi_21 psi_22];
z = zeros(nc+nb,nc+nb);
psi = [psi_1 phi_q';phi_q z];
end