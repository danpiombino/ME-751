function p = iorient(A)
% pi = iorient(A)
% Returns the euler parameters pi for orientation matrix A
%
% Inputs:   A = orientation matrix of reference frame i
%
% Outputs:  p = euler parameter representation of reference frame i
%                   [e0 e1 e2 e3]'
%
% NOTE: DOESN'T WORK FOR 180 DEG ROTATIONS
% Written by: Dan Piombino
% 10/20/16

p(1,1) = sqrt((trace(A)+1)/4);
p(2,1) = (A(3,2)-A(2,3))/4/p(1);
p(3,1) = (A(1,3)-A(3,1))/4/p(1);
p(4,1) = (A(2,1)-A(1,2))/4/p(1);
end