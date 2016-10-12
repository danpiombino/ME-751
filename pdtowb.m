function wi_b = pdtowb(q,q_d,i)
% wi_b = pdtowb(q,q_d,i)
% Returns B matrix for body i given q and ai_b
% 
% Inputs:   q = [7*nb x 1] vector of positions and euler parameters
%               in the form of [r1 r2 ... rnb p1 p2 ... pnb]'
%           i = index of coordinate system of interest
%           q_d = [7*nb x 1] vector of velocities in r and p directions
%
% Outputs:  Ai_d = Time derivative of orientation matrix of coordinate
%                  system i
%
% Written by: Dan Piombino
% 10/12/16

nb = length(q)/7;

if i == 0
    e = [1 0 0 0]'; %by definition
    e_d = zeros(4,1); %by definition
else
    e = q([4*(i-1)+3*nb+1:4*(i-1)+3*nb+4]);
    e_d = q_d([4*(i-1)+3*nb+1:4*(i-1)+3*nb+4]);
end

e_til = skew3(e(2:end));
G = [-e(2:end) -e_til+e(1)*eye(3)];
wi_b = 2*G*e_d;