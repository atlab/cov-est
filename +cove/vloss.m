function L = vloss(R,Rp,V,N)
% normal loss function with condition-specific variances
%
% R :    correlation matrix of training set
% Rp :   correlation matrix of testing set based on training variances
% V :    nBins*nConds*1*nCells - variances per bin per cell
% N :    nBins*nConds*1*nCells - numbers of trials for each element of V

assert(all(size(N)==size(V)))
p = size(R,1);
assert(size(V,4)==p)
N = N/sum(N(:));
L = (trace(Rp/R)+cove.logDet(R))/p + sum(log(V(:)).*N(:));