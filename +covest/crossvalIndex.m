function out = crossvalIndex(n,kfold,chunkSize)
% This function is similar to Matlab's crossvaldin with the 'kfold' option
% except the indices are split into contiguous blocks of approximately
% between chunkSize and 2*chunkSize in length.

assert(chunkSize*kfold*2 <= n, 'chunkSize is too large')
maxBlockSize = ceil(n/kfold);
unassigned = ones(kfold,1)*maxBlockSize;  % unassigned samples for each partition
out = nan(1,n);

fold = randi(kfold);
span = randi(chunkSize);
curr = 0;
while any(unassigned)
    unassigned(fold) = unassigned(fold) - span;
    out(min(curr+(1:span),end)) = fold;
    curr = curr+span;
    
    fold = mod(fold+randi(kfold-1)-1,kfold)+1;
    if unassigned(fold) <= 2*chunkSize
        span = unassigned(fold);
    else
        span = min(ceil(unassigned(fold)/2), chunkSize + randi(chunkSize+1)-1);
    end
end