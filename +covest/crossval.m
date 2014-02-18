function [losses, extras] = crossval(X, loss, estimator, varargin)
% K-fold cross-validation.
% Optimize the hyperparameters, if any, by nested cross validation.
% Return all loss values.
% Compute extras (correlation structure) by applying averaged hyperparameters
% to entire dataset.

K = 10;  % K-fold cross-validation
[n,p] = size(X);
cvInd = covest.crossvalIndex(n,K,round(n/(3*K)));

losses = nan(K,1);
hypers = cell(K,1);
for k=1:K
    XTrain = X(cvInd~=k,:);
    CTest  = covest.cov(X(cvInd==k,:));
    if isempty(varargin)
        C = covest.estimate(XTrain, estimator);
    else
        [C, hypers{k}] = ...
            covest.crossEstimateHyper(XTrain,loss,estimator,varargin{:});
    end
    losses(k) = loss(C,CTest);
end

if nargout>1  % Asking for extras
    if isempty(varargin)
        hypers = [];
    else
        hypers = mean(vertcat(hypers{:}));
    end
    [~,extras] = covest.estimate(X,estimator,hypers);
end