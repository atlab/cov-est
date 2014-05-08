function degree = nodeDegree(S)
% return the average node degree of connectivity in matrix S

p = size(S,1);
[i,j] = meshgrid(1:p,1:p);
degree = sum(abs(S(i~=j))>0)/p;