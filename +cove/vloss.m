function L = vloss(XTest, R, M, V, delta)
% normal loss function based on condition-specific variances
%
% XTest  testing dataset
% R      correlation matrix of training set
% Rp     correlation matrix of testing set based on training variances
% V      nBins*nConds*1*nCells - variances per bin per cell
% N      nBins*nConds*1*nCells - numbers of trials for each element of V in XTest

p = size(XTest,4);
assert(size(V,4)==p)
N = sum(~isnan(XTest),3);
N = N/sum(N(:));
assert(all(size(N)==size(V)))

% subtract means
XTest = bsxfun(@minus, XTest, M);

% shrink the variance estimates toward the mean variance across all conditions
V0 = reshape(nanmean(reshape(V,[],p)),[1 1 1 p]);
V = bsxfun(@plus,(1-delta)*V,delta*V0);

% normalize by regularized variance
Z = bsxfun(@rdivide, XTest, sqrt(V));

% covariance of z-score, i.e. the correlation estimate
Rp = cove.cov(reshape(Z,[],p));

% normal loss
L = (trace(Rp/R)+cove.logDet(R))/p + sum(log(V(:)).*N(:));
if isnan(L)
    L = inf;
end