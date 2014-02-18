function s = sparsity(C)
% return the fraction of zeros in off-diagonal elements of square matric C

p = size(C,1);
[i,j] = meshgrid(1:p,1:p);
s = 1-sum(abs(C(i~=j))>0)/(p*p-p);