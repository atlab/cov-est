function [C, V, extras] = estimate(X, M, varEstimation, corrEstimation, hypers)
% estimate covariance matrix C
%
% Input: 
%    X := nBins * nConds * nTrials * nCells
%    X is zero-mean signal
% 
% Output:
%    C  - correlation matrix
%    V  - matrix of variances:  nBins * nConds * 1 * nCells
%    extras - structure with additional information about the estimate
  

extras = struct;
[nBins, nConds, nTrials, nCells] = size(X);

switch varEstimation
    case 'uniform'
        V = reshape(nanvar(reshape(X,[],nCells)), 1,1,1,nCells);        
        
    case 'linear to mean'
        % variance is assumed to be proportional to mean
        V = nanvar(X,[],3);
        b = arrayfun(@(i) regress(reshape(V(:,:,:,i),[],1),reshape(M(:,:,:,i),[],1)), 1:nCells);
        V = bsxfun(@times, M, reshape(b,1,1,1,nCells));
    
    case 'per condition'
        V = reshape(nanvar(reshape(permute(X, [1 3 2 4]),[],nConds,nCells)), 1,nConds,1,nCells);
        
    case 'per bin'
        V = nanvar(X,[],3);
        
    otherwise 
        error 'unknown variance estimation'    
end

% select cells that have 
    
switch corrEstimation
    case 'sample'
        assert(nargin<3 || isempty(hypers), 'invalid hyperparameters')
        v = sqrt(max(V,0.001*median(V(:))));
        X = bsxfun(@rdivide, X, v);
        C = covest.cov(reshape(X,[],nCells));
        C = C - diag(diag(C)) + eye(size(C));  % replace diagonal with 1s
        
     case 'diag'
         % takes two hypers:  variance shrinkage,  correlation shrinkage
         assert(length(hypers)==2, 'invalid hyperparameters')
         C = covest.shrink(X,hypers(1),hypers(2));
         
%     case 'factor'
%         % takes to hypers: nlatent and shrinkage intensity
%         assert(length(hypers)==2, 'invalid hyperparameters')
%         C = covest.cov(X);
%         [L,psi] = covest.factor(C,hypers(1));
%         C = (1-hypers(2))*C + hypers(2)*(L*L' + diag(psi));  % shrink toward factor model
%         extras.loading_matrix = L;
%         extras.indep_vars = psi;
%         
%     case 'lv-glasso'
%         % set lv-glasso options
%         C = covest.cov(X);
%         scale = sqrt(mean(diag(C)));  % normalize matrix
%         C = scale\C/scale;
%         alpha = hypers(2);
%         beta = hypers(1);
% %        beta = 0;
% %        extras = covest.lvglasso(C, alpha, beta, struct('refit',true,'max_latent',hypers(1)));
%         extras = covest.lvglasso(C, alpha, beta, struct('refit',true,'max_latent',inf));
%         extras.S = scale\extras.S/scale;   % scale back 
%         extras.L = scale\extras.L/scale;
%         [H,D] = svds(extras.L,sum(~~extras.eigL));
%         extras.H = H*sqrt(D);
%         C = inv(extras.S - extras.H*extras.H');
        
    otherwise
        error 'unknown correlation estimator'
end