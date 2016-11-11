function K = kmat(a_b,b)
% Calculates K matrix
% Written by: Dan Piombino
% 11/10/16

K = 2*[a_b'*b a_b'*skew3(b);skew3(a_b)*b a_b*b'+b*a_b'-a_b'*b*eye(3)];
end