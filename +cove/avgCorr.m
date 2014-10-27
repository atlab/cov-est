function m = avgCorr(C)
% mean correlation coefficient
[rows,cols] = size(C);
assert(rows==cols);
%C = corrcov(C);
[i,j] = ndgrid(1:rows,1:cols);
m = mean(C(i<j));
