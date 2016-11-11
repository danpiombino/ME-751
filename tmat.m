function T = tmat(a);
% Calculates T matrix of vector a
% Written by: Dan Piombino
% 11/10/16

T = [0 -a';a -skew3(a)];
end