% log of determinant for large matrices.
% This computation is performed in the log space to avoid num

function [ld,A] = logDet(A)
d = eig(A);
if ~all(d>=0)
    ld = nan;
else
    ld = sum(log(d));
end