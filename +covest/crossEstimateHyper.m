function [C, hypers, extras] = crossEstimateHyper(X, loss, estimator, varargin)
% find the values of hyperparamters that minimize the cross-validated loss
% by K-fold cross-validation.
%
% INPUTS:
%     X - n*p dataset, n=sample size, p=dimension.
%     loss - loss function to be minimized, loss(C,Sigma) where C is estimate
%     estimator -  the name of the estimator as accepted by covest.estimate
%     varargin - lists of valid values for each hyperparameter
%
% OUTPUTS:
%     C - optimized covariance estimate
%     varargout - optimal values of hyperparameters

[n,p] = size(X);
searchSpace = varargin;             % the list of values for each hyperparameter
dims = cellfun(@length, searchSpace);
nHypers = length(dims);
visited = [];  % visited indices
losses =  [];                 % measured losses at visited indices
K = 10;  % K-fold cross-validation
cvInd = covest.crossvalIndex(n,K,round(n/(3*K)));


% coarse random search to seed the local search
disp 'random search'
decimate = 10;
for i=1:prod(ceil(dims/decimate))  % effectively decimate by a factor in each dimension
    visit(arrayfun(@(d) randi(d), dims))
end

% pattern search for the optimum hyperparameter values
disp 'pattern search'
pattern = dims>1;
step = min(floor(dims/2),ceil(dims/decimate*2));
bestLoss = min(losses);
while true
    lastBestLoss = bestLoss;
    [bestLoss,j] = min(losses);
    if bestLoss == lastBestLoss && all(step<=1)
        break
    end
    step = ceil(step/2);
    [best{1:nHypers}] = ind2sub(dims,visited(j));
    ix = [best{:}];
    for b = 0:2^sum(pattern)-1
        s = step;
        s(pattern) = s(pattern).*(2*bitget(b,1:sum(pattern))-1);
        visit(ix+s)
    end
end

[~,j] = min(losses);
[mix{1:nHypers}] = ind2sub(dims,visited(j));
hypers = cellfun(@(h,ix) h(ix), searchSpace, mix);
fprintf('final hyperparameters: %s\n', sprintf(' %g', hypers))
[C,extras] = covest.estimate(X,estimator,hypers);



    function visit(ix)
        % visit a point on the search space specified by ix
        % but return if already visited
        ix = num2cell(max(1,min(dims,ix)));
        hypers_ = cellfun(@(x,i) x(i), searchSpace, ix);
        ix = sub2ind(dims,ix{:});
        alreadyVisited = ismember(ix,visited);
        if ~alreadyVisited
            visited(end+1) = ix;
            losses(end+1) = getLoss(hypers_);
        end
    end


    function L = getLoss(hypers)
        fprintf('%g  ', hypers)
        L = nan(1,K);
        for k=1:K
            % record average loss
            CC = covest.estimate(X(cvInd~=k,:),estimator, hypers);
            testC = covest.cov(X(cvInd==k,:));
            L(k) = loss(CC,testC);
        end
        if any(isnan(L))
            L = inf;
        else
            L = mean(L);
        end
        fprintf(':  mean loss %g\n', L)
    end
end