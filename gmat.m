function G = gmat(q,i)
% G = gmat(q,i)
% Returns G matrix for body i given q
% 
% Inputs:   q = [7*nb x 1] vector of positions and euler parameters
%               in the form of [r1 r2 ... rnb p1 p2 ... pnb]'
%           i = index of coordinate system of interest (0 = ground)
%
% Outputs:  G = G matrix of coordinate system i
%
% Note:     Input q_d instead of q to obtain G_d
%
% Written by: Dan Piombino
% 10/28/16

nb = length(q)/7;

if i == 0
    e = [1 0 0 0]'; %by definition
else
    e = q([4*(i-1)+3*nb+1:4*(i-1)+3*nb+4]);
end
G = [-e(2:4), -skew3(e(2:4))+e(1)*eye(3)];
end