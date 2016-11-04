function [r,p] = q2rp(q)
% rearranges from [r1 p1 r2 p2 ...]' to [r1 r2 r3... p1 p2 p3...]
% Written by Dan Piombino
% 11/3/16

nb = length(q)/7;
r = zeros(3*nb,1);
p = zeros(4*nb,1);
for j = 1:nb
    r(3*(j-1)+1:3*(j-1)+3,1) = q(7*(j-1)+1:7*(j-1)+3,1);
    p(4*(j-1)+1:4*(j-1)+4,1) = q(7*(j-1)+4:7*j,1);
end
end