function [C, extras] = estimate(X,estimator,hypers)
% estimate covariance matrix C
% some estimators return additional results in extras.
extras = struct;

switch estimator
    case 'sample'
        assert(nargin<3 || isempty(hypers), 'invalid hyperparameters')
        C = covest.cov(X);
        
    case 'shrink'
        % takes two hypers:  variance shrinkage,  correlation shrinkage
        assert(length(hypers)==2, 'invalid hyperparameters')
        C = covest.shrink(X,hypers(1),hypers(2));
        
    case 'factor'
        % takes to hypers: nlatent and shrinkage intensity
        assert(length(hypers)==2, 'invalid hyperparameters')
        C = covest.cov(X);
        [L,psi] = covest.factor(C,hypers(1));
        C = (1-hypers(2))*C + hypers(2)*(L*L' + diag(psi));  % shrink toward factor model
        extras.loading_matrix = L;
        extras.indep_vars = psi;
        
    case 'lv-glasso'
        % set lv-glasso options
        C = covest.cov(X);
        scale = sqrt(mean(diag(C)));  % normalize matrix
        C = scale\C/scale;
        alpha = hypers(2);
        beta = hypers(1);
%        beta = 0;
%        extras = covest.lvglasso(C, alpha, beta, struct('refit',true,'max_latent',hypers(1)));
        extras = covest.lvglasso(C, alpha, beta, struct('refit',true,'max_latent',inf));
        extras.S = scale\extras.S/scale;   % scale back 
        extras.L = scale\extras.L/scale;
        [H,D] = svds(extras.L,sum(~~extras.eigL));
        extras.H = H*sqrt(D);
        C = inv(extras.S - extras.H*extras.H');
        
    otherwise
        error 'unknown estimator'
end