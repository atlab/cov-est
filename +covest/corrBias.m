function corrBias
% assess effect of the common variance assumption on limited range
% correlations

C1 = zeros(0,3);  % average correlation, binwise, condition-wise, and common
C2 = zeros(0,3);  % average correlation of similarly tuned cells
oriRel = pro(covest.OriTuning,'high_repeats->temp','*');
for key = fetch(covest.ActiveCells*oriRel & 'ncells>100')'
    [X,evokedBins,sel,ori,pval] = fetch1(covest.ActiveCells*covest.Traces*oriRel & key, ...
        'trace_segments','evoked_bins', 'selection','von_pref','von_p_value');
    X = double(X(1:min(end,evokedBins),:,:,sel));  % only evoked bins
    ori = ori(sel)*180/pi;
    pval = pval(sel);
    c0 = computeCorr(X,'bin');
    c1 = computeCorr(X,'cond');
    c2 = computeCorr(X,'common');
    
    
    q = [nan 1];
    ori(pval>0.05)=nan;
    od = oriDiff(ori,ori');

    % average correlations between differently tuned cells
    ix = q((~isnan(od) & od > 15)+1);
    [C1(end+1,:),n] = cellfun(@(c) avgOffDiag(c.*ix), {c0 c1 c2});    
    
    % average correlations between similarly tuned cells
    ix = q((~isnan(od) & od < 15)+1);
    [C2(end+1,:),n] = cellfun(@(c) avgOffDiag(c.*ix), {c0 c1 c2});
    
end
fig = Figure(1,'size',[40 40])
bar(mean(C2-C1)')
colormap(0.5+0.5*gray)
hold on
errorbar(mean(C2-C1)',sqrt(var(C1-C2)/size(C1,1)),'k.')
end


function C = corrVar(Z,varComp)
% compute the correlation matrix with different variance assumptions
[nBins,nConds,nTrials,nCells] = size(Z);

M = nanmean(Z,3);
Z = bsxfun(@minus, Z, M);

% common variance
V0 = nanvar(reshape(Z,[],nCells));
V0 = reshape(V0,1,1,1,nCells);
a = 0.01; % very modest regularization

switch varComp
    case 'common'
        V = V0;
    case 'bin'
        % bin-wise variance
        V = nanvar(Z,[],3);
    case 'cond'
        V = mean(nanvar(Z,[],3));
    otherwise
        error 'unknown variance computation'
end
V = bsxfun(@plus, a*V0, (1-a)*V);
Z = bsxfun(@rdivide, Z, sqrt(V));
C = corrcov(cove.cov(reshape(Z,[],nCells)));
end


function d=oriDiff(ori1,ori2)
% compute the absolute difference between orientations ori1 and ori2 (in degrees)
ori1 = mod(ori1, 180);
ori2 = mod(ori2, 180);
s = bsxfun(@plus, ori1, ori2); % ensure that nans remain in place
b1 = s - bsxfun(@max,ori1,ori2);
b2 = s - bsxfun(@min,ori1,ori2);
d = min(b2-b1,b1+180-b2);
end


function [m,n] = avgOffDiag(C)
% mean correlation coefficient
p = size(C,1);
[i,j] = ndgrid(1:p,1:p);
m = nanmean(C(i<j));
n = sum(~isnan(C(i<j)));
end