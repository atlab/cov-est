function fig6_ori
close all

% get all data
c = covest.CovMatrix & 'nfolds=1' & 'sparsity<0.97';
c0 = c.pro('method->m0', 'corr_matrix->c0') & 'm0=0';
c1 = c.pro('method->m1','sparse','corr_matrix','sparsity') & 'm1=91';
rel = c0*c1*covest.ActiveCells*pro(covest.OriTuning,'high_repeats->hrepeats','*');
[pval,ori,selection,C0,C1,S,sparsity] = ...
    rel.fetchn('von_p_value','von_pref','selection','c0','corr_matrix','sparse','sparsity');

% remove inactive cells
ori = cellfun(@(ori,selection) ori(selection)*180/pi, ori, selection,'uni',false);
pval = cellfun(@(pval,selection) pval(selection), pval, selection,'uni',false);
clear selection

% convert to correlations and partial correlations, respectively
C0 = cellfun(@corrcov, C0, 'uni', false);  % sample correlations
C1 = cellfun(@(C) -corrcov(inv(C)), C1, 'uni', false);  % regularized partial correlations

% remove untuned cells
alpha = .05;
ori = cellfun(@(ori,pval) ori(pval<alpha), ori, pval, 'uni',false);
C0 = cellfun(@(C,pval) C(pval<alpha,pval<alpha), C0, pval, 'uni',false);
C1 = cellfun(@(C,pval) C(pval<alpha,pval<alpha), C1, pval, 'uni',false);
S  = cellfun(@(C,pval) C(pval<alpha,pval<alpha), S, pval, 'uni',false);
clear pval

% matrices of orientation differences
D = cellfun(@oriDiffMatrix, ori, 'uni',false);

% get off-diagonal elements
D  = cellfun(@offDiag, D,  'uni', false);
C0 = cellfun(@offDiag, C0, 'uni',false);
C1 = cellfun(@offDiag, C1, 'uni',false);
S  = cellfun(@offDiag, S,  'uni',false);

% bin orientation differences
edges = [0 15 45 90];
nBins = length(edges)-1;
D = cellfun(@(oriDiffs) sum(bsxfun(@ge,oriDiffs,edges),2), D, 'uni', false);
 
doStats = true;
% do statistics
if doStats
    ds = cellfun(@(D,S) mat2dataset([D S>0 S<0], 'VarNames', {'distance','neg','pos'}), D, S, 'uni', false);
    mdlPos = cellfun(@(ds) fitglm(ds,'pos ~ distance','Distribution','binomial'), ds, 'uni', false);
    mdlNeg = cellfun(@(ds) fitglm(ds,'neg ~ distance','Distribution','binomial'), ds, 'uni', false);
    posP = cellfun(@(mdl) mdl.Coefficients.pValue(2), mdlPos);
    negP = cellfun(@(mdl) mdl.Coefficients.pValue(2), mdlNeg);
    
    tstat = cellfun(@(pos,neg) (pos.Coefficients.Estimate(2)-neg.Coefficients.Estimate(2))/norm([pos.Coefficients.SE(2) neg.Coefficients.SE(2)]), mdlPos, mdlNeg);
    diffP = arrayfun(@(t,n) tcdf(t,n), tstat, cellfun(@(D) sum(D>0),D));  % p-value of difference between the coefficients
end

% exclude sites that don't have enough tuned cell pairs in each bin
hasEnough = 200 < cellfun(@(D) min(hist(D,1:nBins)), D);
fprintf('Qualifying sites n=%d\n', sum(hasEnough))
D = D(hasEnough);
C0 = C0(hasEnough);
C1 = C1(hasEnough);
S = S(hasEnough);
sparsity = sparsity(hasEnough);
clear hasEnough

% bin correlations into bins specified by edges
C0 = cellfun(@(C,D) accumarray(D,C,[nBins 1],@mean), C0, D, 'uni',false);
C1 = cellfun(@(C,D) accumarray(D,C,[nBins 1],@mean), C1, D, 'uni',false);
C0 = cat(2,C0{:})';
C1 = cat(2,C1{:})';

% panel A: average correlations vs ori tuning
fig = Figure(1,'size',[45 32]);
ticks = conv(edges,[.5 .5],'valid')';
plot(ticks,C0,'k.-','LineWidth',.25);
yticks = 0:0.05:1;
set(gca,'XTick',edges,'YTick',yticks,'YTickLabel',nozero(yticks))
xlim([0 90])
ylim([0 0.055])
xlabel '\Deltaori (degrees)'
ylabel 'avg. noise corr.'
set(gca,'Position',[0.27 0.29 0.67 0.67])
fig.cleanup
fig.save(fullfile(covest.plots.figPath, 'Fig6-A.eps'))

% panel D: average partial correlations vs ori tuning
fig = Figure(1,'size',[45 32]);
plot(ticks,C1,'k.-','Linewidth',.25);
xlim([0 90])
ylim([0 0.014])
yticks = 0:0.01:0.02;
set(gca,'XTick',edges,'YTick',yticks,'YTickLabel',arrayfun(@(s) strrep(sprintf('%g',s),'0.','.'), yticks,'uni',false))
xlabel '\Deltaori (degrees)'
ylabel 'avg. partial corr.'
set(gca,'Position',[0.27 0.29 0.67 0.67])
fig.cleanup
fig.save(fullfile(covest.plots.figPath, 'Fig6-D.eps'))

% panel G: connectivity vs ori tuning
fig = Figure(1,'size',[45 32]);
posColor = [0 .5 0];
negColor = [.5 0 0];

neg = cellfun(@(S,D) accumarray(D,S>0,[nBins 1], @mean), S, D, 'uni',false);
pos = cellfun(@(S,D) accumarray(D,S<0,[nBins 1], @mean), S, D, 'uni',false);

plot(ticks,cat(2,pos{:}),'^-', 'Color',posColor,'MarkerSize',3,'LineWidth',.25)
hold on
plot(ticks,cat(2,neg{:}),'v:', 'Color',negColor,'MarkerSize',3,'LineWidth',.25)
hold off
yticks = 0:0.1:0.2;
set(gca,'XTick',edges,'YTick',yticks,'YTickLabel',arrayfun(@(s) strrep(sprintf('%g',s),'0.','.'), yticks,'uni',false))
xlim([0 90])
ylim([0 0.25])
xlabel '\Deltaori (degrees)'
ylabel 'connectivity'
set(gca,'Position',[0.27 0.29 0.67 0.67])
fig.cleanup
fig.save(fullfile(covest.plots.figPath, 'Fig6-G.eps'))
end


function ret = oriDiffMatrix(ori)
% make a matrix of orientation differences
ret = oriDiff(ori,ori');
end


function d=oriDiff(ori1,ori2)
% compute the absolute difference between orientations ori1 and ori2 (in degrees)
ori1 = mod(ori1, 180);
ori2 = mod(ori2, 180);
s = bsxfun(@plus, ori1, ori2); % ensure that nans remain in place since max(nan,A)==A
b1 = s - bsxfun(@max,ori1,ori2);
b2 = s - bsxfun(@min,ori1,ori2);
d = min(b2-b1,b1+180-b2);
end

function ret = offDiag(C)
% return the vector of off-diagonal elements of symmetric matrix C
p = size(C,1);
[i,j] = ndgrid(1:p,1:p);
ret = C(i<j);
end


function s=nozero(f)
% remove leading zeros in decimal fractions
s = arrayfun(@(f) strrep(sprintf('%g',f),'0.','.'), f, 'uni', false);
end