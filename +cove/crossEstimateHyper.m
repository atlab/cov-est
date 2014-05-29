function [hypers, bestDelta, visited, losses] = crossEstimateHyper(X, evokedBins, reg, searchSpace)
% find the values of hyperparamters that minimize the cross-validated loss
% by K-fold cross-validation.
%
% INPUTS:
%     X = nBins * nDirs * nTrials * nCells
%     covEstimation -  string identifying the method for estimating the variance
%     reg - structure specifying regularization of means, variances, and correlations
%     searchSpace - lists of valid values for each hyperparameter, in sequence
%
% OUTPUTS:
%     hypers - optimal values of hyperparameters
%     visited - the indices of visited hyperparameter values
%     losses  - the loss function values for these hyperparameter values

dims = cellfun(@length, searchSpace);
nHypers = length(dims);
assert(nHypers>0)

visited = [];  % visited indices
losses =  [];  % measured losses at visited indices

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
    assert(~isnan(bestLoss))
    if bestLoss == lastBestLoss && all(step<=1)
        break
    end
    step = ceil(step/2);
    [best{1:nHypers}] = ind2sub(dims,visited(j));
    ix = [best{:}];
    % visit all nodes 1 step away in every direction
    for b = 0:2^sum(pattern)-1
        s = step;
        s(pattern) = s(pattern).*(2*bitget(b,1:sum(pattern))-1);
        visit(ix+s)
    end
end
[~,j] = min(losses);
[indices{1:nHypers}] = ind2sub(dims,visited(j));
hypers = cellfun(@(h,ix) h(ix), searchSpace, indices);
fprintf('final hyperparameters: %s\n', sprintf(' %g', hypers))


    function visit(ix)
        % visit a point on the search space specified by ix
        % but return if already visited

        K = 10;        % K-fold cross-validation

        ix = num2cell(max(1,min(dims,ix)));
        hypers_ = cellfun(@(x,i) x(i), searchSpace, ix);
        ix = sub2ind([dims ones(1,2-length(dims))],ix{:});
        alreadyVisited = ismember(ix,visited);
        if ~alreadyVisited
            fprintf('%2.4g  ', hypers_)
            [XTest,R,M,V] = arrayfun(@(k) estimate_(hypers_,k,K), 1:K, 'uni', false);
            delta = mean(cellfun(@(XTest,R,M,V) cove.findBestDelta(XTest, R, M, V), XTest, R, M, V));
            visited(end+1) = ix; 
            losses(end+1) = mean(cellfun(@(XTest,R,M,V) cove.vloss(XTest,R,M,V,delta), XTest, R, M, V));
            if losses(end)==min(losses)
                bestDelta = delta;
            end
            fprintf(':  mean loss %g\n', losses(end))
        end
    end

    function [XTest,R,M,V] = estimate_(hypers,k,K)
        % compute cross-validation loss
        [XTrain,XTest] = cove.splitTrials(X,k,K);
        [R, M, V] = cove.estimate(XTrain, evokedBins, reg, hypers);        
    end
end