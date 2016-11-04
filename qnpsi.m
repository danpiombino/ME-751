function psi = qnpsi(M,Jp,phi_q)
% function psi = qnpsi(M,Jp,phi_q)
% Calculates psi for a quasi-newton method
% 
% Inputs:   M: Mass matrix [3*nb x 3*nb]
%           Jp: Euler Parameter Inertia Matrix [4*nb x 4*nb]
%           phi_q: Jacobian of constraint functions [nc+nb x 7*nb]
%
% Outputs:   psi: Jacobian of g(z) for quasi-newton method
% 
% Written By: Dan Piombino
% 11/3/16

nb = size(M,1)/3;
nc = size(phi_q,1)-nb;

z12 = zeros(3*nb,4*nb);
z14 = zeros(3*nb,nb);
z21 = zeros(4*nb,3*nb);
z33 = zeros(nc,nc);
z34 = zeros(nc,nb);
z41 = zeros(nb,3*nb);
z43 = zeros(nb,nc);
z44 = zeros(nb,nb);

phi_r = phi_q(1:nc,1:3*nb);
phi_p = phi_q(1:nc,3*nb+1:end);
P = phi_q(nc+1:end,3*nb+1:end);

psi = [M z12 phi_r' z14;z21 Jp phi_p' P';phi_r phi_p z33 z34;z41 P z43 z44];
end