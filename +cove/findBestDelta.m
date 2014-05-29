function delta = findBestDelta(XTest, R, M, V)
% search for optimal variance regularization toward common variance
% across all conditions
delta = exp(-4:.05:0);
L = nan(size(delta));
idx = floor(length(delta)/2)+(1:2);
for i=idx
    L(i) = cove.vloss(XTest, R, M, V, delta(i));
end
assert(~any(isinf(L(idx))))
if L(idx(1))<L(idx(2))
    idx = idx(1)-1:-1:1;
else
    idx = idx(2)+1:length(delta);
end
for i=idx
    L(i) = cove.vloss(XTest, R, M, V, delta(i));
    if L(i)>min(L)
        break
    end
end
delta = delta(i);