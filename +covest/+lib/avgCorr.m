function m = avgCorr(C)
p = size(C,1);
C = corrcov(C);
[i,j] = ndgrid(1:p,1:p);
m = mean(C(i<j));