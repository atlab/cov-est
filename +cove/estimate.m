function [R, M, V, extras] = estimate(X, evokedBins, reg, hypers)
% estimate covariance matrix C
%
% Input:
%    X - data: nBins * nConds * nTrials * nCells
%    evokedBins - the number of bins considered in the evoked response,
%    whereas the remaning bins are considered to be spontaneous
%    reg - covariance regularization method
%
% Output:
%    R  - correlation matrix
%    M  - mean responses:  nBins * nConds * nTrials * nCells
%    V  - variances, dimensions:  nBins * nConds * 1 * nCells
%    extras - structure with additional information about the estimate,
%    i.e. sparse component, low-rank components, sparsity, etc.

extras = struct;
[nBins, nConds, nTrials, nCells] = size(X);
evokedBins = min(size(X,3),evokedBins);

% estimate mean response
% bin-specific means
M = nanmean(X(1:evokedBins,:,:,:),3);   % binwise mean of evoked response
if evokedBins < size(X,1)
    % common mean during spontaneous activity
    M2 = reshape(nanmean(reshape(X(evokedBins+1:end,:,:,:),[],nCells)),1,1,1,nCells);  % common mean in intertrial periods
    M2 = repmat(M2,size(X,1)-evokedBins,nConds);
    M = cat(1,M,M2);
end

% subtract mean
X = bsxfun(@minus, X, M);

% estimate binwise variances
V = nanvar(X(1:evokedBins,:,:,:),1,3);   % binwise mean of evoked response
if evokedBins < size(X,1)
    % common variance during spontaneous activity
    V2 = reshape(nanvar(reshape(X(evokedBins+1:end,:,:,:),[],nCells),1),1,1,1,nCells);  % common mean in intertrial periods
    V2 = repmat(V2,size(X,1)-evokedBins,nConds);
    V = cat(1,V,V2);
end

% sample correlation matrix based on bin-wise variances
Z = bsxfun(@rdivide, X, sqrt(V));
R = cove.cov(reshape(Z,[],nCells));
assert(all(abs(diag(R)-1)<1e-3))
R = corrcov(R);  % just in case

% average variances across all bins and produce the average
% sample covariance matrix, C
sigma = diag(sqrt(mean(reshape(V,[],nCells))));  % average variances
C = sigma*R*sigma;

switch reg
    case 'sample'
        assert(isempty(hypers),'invalid hyperparameters')
        % do nothing
        
    case 'diag'
        % 2 hypers:  variance shrinkage,  correlation shrinkage
        assert(length(hypers)==2, 'invalid hyperparameters')
        v = diag(C);
        v = diag(sqrt((1-hypers(1))*v + hypers(1)*median(v)));
        r = (1-hypers(2))*corrcov(C) + hypers(2)*eye(size(C));
        C = v*r*v;
        
    case 'factor'
        % 2 hypers: variance shrink toward median, nlatent
        assert(length(hypers)==2, 'invalid hyperparameters')
        v = diag(C);
        v = diag(sqrt((1-hypers(1))*v + hypers(1)*median(v)));
        C = v*corrcov(C)*v;
        [L,psi] = cove.factor(C,hypers(2));  % factor analysis
        extras.loading_matrix = L;
        extras.indep_vars = psi(:);
        
    case 'glasso'
        assert(length(hypers)==1)
        cove.set('max_latent',0)   % prevent latent variables
        scale = mean(diag(C));
        extras = cove.lvglasso(C/scale,hypers(1),10,cove.set);
        extras.S = extras.S/scale;
        C = inv(extras.S);
        
    case 'lv-glasso'
        assert(length(hypers)==2)
        cove.set('max_latent',inf)   % prevent latent variables
        scale = mean(diag(C));
        extras = cove.lvglasso(C/scale,hypers(1),hypers(2),cove.set);
        extras.S = extras.S/scale;  % scale back
        extras.L = extras.L/scale;
        [H,D] = svds(extras.L,sum(~~extras.eigL));
        extras.H = H*sqrt(D);
        C = inv(extras.S - extras.H*extras.H');
        
    otherwise
        error 'unknown covariance estimator'
end

% convert back to correlations
R = sigma\C/sigma;

if ~strcmp(reg,'sample')
    % transfer the change in variance from R to V
    V = bsxfun(@times, V, reshape(diag(R), [1 1 1 nCells])); 
    R = corrcov(R);
end