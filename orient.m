function Ai = orient(q,i)
% Ai = orient(q,i)
% Returns the orientation matrix A for body i given q
%
% Inputs:   q = [7*nb x 1] vector of positions and euler parameters
%               in the form of [r1 r2 ... rnb p1 p2 ... pnb]'
%           i = index of coordinate system of interest (0 = ground)
%
% Outputs:  Ai = orientation matrix of coordinate system i with respect to ground
%
% Written by: Dan Piombino
% 10/12/16

if i == 0;
    Ai = eye(3);
else
    nb = length(q)/7;
    e = q([4*(i-1)+3*nb+1:4*(i-1)+3*nb+4]);
    Ai(1,1) = e(1)^2+e(2)^2-0.5;
    Ai(1,2) = e(2)*e(3)-e(1)*e(4);
    Ai(1,3) = e(2)*e(4)+e(1)*e(3);
    Ai(2,1) = e(2)*e(3)+e(1)*e(4);
    Ai(2,2) = e(1)^2+e(3)^2-0.5;
    Ai(2,3) = e(3)*e(4)-e(1)*e(2);
    Ai(3,1) = e(2)*e(4)-e(1)*e(3);
    Ai(3,2) = e(3)*e(4)+e(1)*e(2);
    Ai(3,3) = e(1)^2+e(4)^2-0.5;
    Ai = 2*Ai;
end
end