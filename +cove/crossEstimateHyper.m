function [hypers, visited, losses] = crossEstimateHyper(X, evokedBins, loss, covEstimation, searchSpace)
% find the values of hyperparamters that minimize the cross-validated loss
% by K-fold cross-validation.
%
% INPUTS:
%     X = nBins * nDirs * nTrials * nCells
%     loss - loss function to be minimized, loss(C,Sigma) where C is estimate
%     covEstimation -  string identifying the method for estimating the variance
%     hypers - lists of valid values for each hyperparameter, in sequence
%
% OUTPUTS:
%     C - optimized covariance estimate
%     varargout - optimal values of hyperparameters

dims = cellfun(@length, searchSpace);
nHypers = length(dims);
assert(nHypers>0)

visited = [];  % visited indices
losses =  [];  % measured losses at visited indices
K = 10;        % K-fold cross-validation

% coarse random search to seed the local search
fprintf 'random search: '
decimate = 6;
nRandom = ceil(prod(max(1,dims/decimate)/sum(dims>1)));
fprintf('%d points\n', nRandom)
ix = arrayfun(@(~) arrayfun(@(d) randi(d), dims), 1:nRandom, 'uni', false);
cellfun(@visit,ix);

% pattern search for the optimum hyperparameter values
disp 'pattern search'
pattern = dims>1;
step = min(floor(dims/2),ceil(dims/decimate*2));
bestLoss = min(losses);
while true
    lastBestLoss = bestLoss;
    [bestLoss,j] = min(losses);
    if isnan(bestLoss) || bestLoss == lastBestLoss && all(step<=1)
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



    function visit(ix)
        % visit a point on the search space specified by ix
        % but return if already visited
        ix = num2cell(max(1,min(dims,ix)));
        hypers_ = cellfun(@(x,i) x(i), searchSpace, ix);
        ix = sub2ind([dims ones(1,2-length(dims))],ix{:});
        alreadyVisited = ismember(ix,visited);
        if ~alreadyVisited
            visited(end+1) = ix;
            losses(end+1) = cvLoss(hypers_);
        end
    end


    function L = cvLoss(hypers)
        % compute avearge cross-validation loss
        fprintf('%g  ', hypers)
        L = nan(1,K);
        for k=1:K
            % record average loss
            [XTrain,XTest] = cove.splitTrials(X,k,K);
            C0 = cove.estimate(XTrain, [], evokedBins, 'sample', {});
            C = cove.estimate(XTrain, [], evokedBins, covEstimation, hypers);
            L(k) = loss(C, ...
                cove.estimate(XTest, [], evokedBins, 'sample', {}));
            if isnan(L(k)), break, end
        end
        L = mean(L);
        if isnan(L)
            L = inf;
        end
        fprintf(':  mean loss %g\n', L)
    end
end