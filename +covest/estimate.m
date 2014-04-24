function [C, M, extras] = estimate(X, M, evokedBins, covEstimation, hypers)
% estimate covariance matrix C
%
% Input:
%    X := nBins * nConds * nTrials * nCells
%
% Output:
%    C  - covariance matrix
%    extras - structure with additional information about the estimate


extras = struct;
[nBins, nConds, nTrials, nCells] = size(X);

% subtract mean response
if isempty(M)
    if size(X,1)<=evokedBins
        M = nanmean(X,3);
    else
        M1 = nanmean(X(1:evokedBins,:,:,:),3);   % binwise mean of evoked response
        M2 = reshape(nanmean(reshape(X(evokedBins+1:end,:,:,:),[],nCells)),1,1,1,nCells);  % common mean in intertrial periods
        M2 = repmat(M2,size(X,1)-evokedBins,nConds);
        M = cat(1,M1,M2);
    end
end

X = bsxfun(@minus, X, M);
X = reshape(X,[],nCells);
X = X(~any(isnan(X),2),:);


switch covEstimation
    case 'sample'
        assert(isempty(hypers),'invalid hyperparameters')
        C = covest.lib.cov(X);
        
    case 'diag'
        % takes two hypers:  variance shrinkage,  correlation shrinkage
        assert(length(hypers)==2, 'invalid hyperparameters')
        C = covest.lib.shrink(X,hypers(1),hypers(2));
        
    case 'factor'
        % takes 3 hypers: nlatent, shrink, and individual variance shrink
        assert(length(hypers)==3, 'invalid hyperparameters')
        C = covest.lib.cov(X);
        [L,psi] = covest.lib.factor(C,hypers(1));
        psi = (1-hypers(3))*psi + hypers(3)*mean(psi);
        C = (1-hypers(2))*C + hypers(2)*(L*L' + diag(psi));  % shrink toward factor model
        extras.loading_matrix = L;
        extras.indep_vars = psi;
        
    case 'glasso'
        assert(length(hypers)==2)
        C = covest.lib.cov(X);
        C = hypers(1)*mean(diag(C))*eye(size(C)) + (1-hypers(1))*C; 
        
        scale = mean(diag(C));
        C = C/scale;
        alpha = hypers(2);
        beta = 10;
        extras = covest.lib.lvglasso(C,alpha,beta,struct('refit',true,'max_latent',0));
        extras.S = extras.S/scale;  % scale back
        C = inv(extras.S);
        
        
    case 'lv-glasso'
        % set lv-glasso options
        C = covest.lib.cov(X);
        C = hypers(1)*mean(diag(C))*eye(size(C)) + (1-hypers(1))*C;
        
        scale = mean(diag(C));
        C = C/scale;
        alpha = hypers(2);
        beta = hypers(3);
        extras = covest.lib.lvglasso(C,alpha,beta,struct('refit',true));
        extras.S = extras.S/scale;  % scale back
        extras.L = extras.L/scale;
        [H,D] = svds(extras.L,sum(~~extras.eigL));
        extras.H = H*sqrt(D);
        C = inv(extras.S - extras.H*extras.H');
        
    otherwise
        error 'unknown covariance estimator'
end