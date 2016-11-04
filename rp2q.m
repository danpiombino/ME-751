function q = rp2q(r,p)
% rearranges from [r1 r2 r3... p1 p2 p3...] to [r1 p1 r2 p2 ...]'
% Written by Dan Piombino
% 11/3/16

nb = length(r)/3;
q = zeros(7*nb,1);
for j = 1:nb
    q(7*(j-1)+1:7*(j-1)+3,1) = r(3*(j-1)+1:3*(j-1)+3,1);
    q(7*(j-1)+4:7*j,1) = p(4*(j-1)+1:4*(j-1)+4,1);
end
end