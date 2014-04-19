function [X, XTest] = splitTrials(X, fold, nFolds)
% split X's trials into training dataset and testing dataset
% The dimensions of X are: nBins * nConds * nTrials * nCells

nTrials = size(X,3);
if nFolds == 1
    XTest = [];
else
    iTrials = mod(1:nTrials,nFolds)+1;
    rng(0)   % ensure the same folds every time
    iTrials = iTrials(randperm(end));
    XTest = X(:,:,iTrials==fold,:);
    X = X(:,:,iTrials~=fold,:);
end