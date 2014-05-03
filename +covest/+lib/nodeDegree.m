function degree = nodeDegree(C)
% return the average node degree of connectivity in matrix C

p = size(C,1);
[i,j] = meshgrid(1:p,1:p);
degree = sum(abs(C(i~=j))>0)/p;