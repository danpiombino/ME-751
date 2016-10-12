function atil = skew3(a)
% atil = skew3(a)
% Returns the skew-symmetric matrix atil created by the generator vector a
%
% Inputs:   a [column or row vector of length 3]
%
% Outputs:  atil [3 by 3 skew symmetric matrix]
%
% Written by: Dan Piombino
% 10/12/2016

atil = [0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
end