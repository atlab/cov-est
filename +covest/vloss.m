function L = vloss(X,C,V,loss)
% compute the average log likelihood of trials in X given the correlation
% matrix C and variances V. 
%
% X and V dimensions nBins * nConds * nTrials * nCells
% Singleton dimensions in V signify common variance across all values of
% that dimension.

[nBins, nConds, ~, nCells] = size(X);
assert(all(size(V)==[nBins,nConds,1,nCells] | size(V)==[1 1 1 nCells]))
assert(all(size(C)==[nCells nCells]))

% average loss across all conditions, assuming balanced design
L = 0;
totalTrials = 0;
for iBin = 1:size(V,1)
    if size(V,1)==1
        iBin = 1:size(X,1); %#ok<FXSET>
    end
    for iCond = 1:size(V,2)
        if size(V,1)==1
            iCond = 1:size(X,2); %#ok<FXSET>
        end
        v = sqrt(squeeze(V(iBin(1),iCond(1),:,:)));
        v = diag(max(0.001*median(v),v));
        assert(length(v)==nCells);
        [c,nTrials] = covest.cov(reshape(X(iBin,iCond,:,:),[],nCells));
        L = L + loss(v*C*v, c)*nTrials;
        totalTrials = totalTrials + nTrials;
    end
end
L = L / totalTrials;