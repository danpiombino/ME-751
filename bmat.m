function Bi = bmat(q,i,ai_b)
% Bi = bmat(q,i,ai_b)
% Returns B matrix for body i given q and ai_b
% 
% Inputs:   q = [7*nb x 1] vector of positions and euler parameters
%               in the form of [r1 r2 ... rnb p1 p2 ... pnb]'
%           i = index of coordinate system of interest (0 = ground)
%           ai_b = vector a in the local coordinate system i
%
% Outputs:  Bi = B matrix of coordinate system i
%
% EQ used:  B(p,ai_b) = 2*[(e0*I3+e_til)*ai_b,
% e*ai_b'-(e0*I3+e_til)*ai_b_til]
%
% Note:     e = [e1 e2 e3]' and e_til is skew3(e) in the above EQ
%
% Written by: Dan Piombino
% 10/12/16

nb = length(q)/7;

if i == 0
    e = [1 0 0 0]'; %by definition
else
    e = q([4*(i-1)+3*nb+1:4*(i-1)+3*nb+4]);
end
nb = length(q)/7;
e_til = skew3(e(2:end));
ai_b_til = skew3(ai_b);
Bi(1,1) = 2*(e(1)*eye(3)+e_til)*ai_b;
Bi(1,2) = 2*(e(2:end)*ai_b'-(e(1)*eye(3)+e_til)*ai_b_til);
end