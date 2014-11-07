function fig6_dist
close all
% Figure 6, panels B, C, E, F
for vertical = [false true]  % false = panels B, D,  false = C, F
    c = covest.CovMatrix & 'nfolds=1' & 'sparsity<0.97';
    c0 = c.pro('method->m0','corr_matrix->c0') & 'm0=0';
    c1 = c.pro('method->m1','corr_matrix','sparse','sparsity') & 'm1=91';
    rel = c0*c1*covest.Traces*covest.ActiveCells;
    [xyz,selection,C0,C1,S,sparsity] = ...
        rel.fetchn('cell_xyz','selection','c0','corr_matrix','sparse','sparsity');
    
    % remove invactive cells
    xyz = cellfun(@(xyz,selection) xyz(selection,:), xyz, selection, 'uni', false);
    
    % convert to correlations and partial correlations, respectively
    C0 = cellfun(@corrcov, C0, 'uni', false);  % sample correlations
    C1 = cellfun(@(C) -corrcov(inv(C)), C1, 'uni', false);  % regularized partial correlations
    
    % matrices of lateral and vertical distances
    Dx = cellfun(@(xyz) distMatrix(xyz(:,1:2)), xyz, 'uni',false);
    Dz = cellfun(@(xyz) distMatrix(xyz(:,  3)), xyz, 'uni',false);
    
    % get off-diagonal elements
    C0 = cellfun(@offDiag, C0, 'uni',false);
    C1 = cellfun(@offDiag, C1, 'uni',false);
    S  = cellfun(@offDiag, S,  'uni',false);
    Dx = cellfun(@offDiag, Dx, 'uni',false);
    Dz = cellfun(@offDiag, Dz, 'uni',false);
    
    % eliminate cell pairs that are too distant in the orthogonal direction
    maxDist = 30; % um
    if vertical
        edges = [0 25 60 130];
        r = Dx;
        D = Dz;
    else
        edges = [0 25 75 150 250];
        r = Dz;
        D = Dx;
    end
    C0 = cellfun(@(C,r) C(r<maxDist), C0, r, 'uni',false);
    C1 = cellfun(@(C,r) C(r<maxDist), C1, r, 'uni',false);
    S  = cellfun(@(C,r) C(r<maxDist), S,  r, 'uni',false);
    D  = cellfun(@(C,r) C(r<maxDist), D,  r, 'uni',false);
    
    % do statistics
    doStats = false;
    if doStats
        ds = cellfun(@(D,S) mat2dataset([D S>0 S<0], 'VarNames', {'distance','neg','pos'}), D, S, 'uni', false);
        mdlPos = cellfun(@(ds) fitglm(ds,'pos ~ distance','Distribution','binomial'), ds, 'uni', false);
        mdlNeg = cellfun(@(ds) fitglm(ds,'neg ~ distance','Distribution','binomial'), ds, 'uni', false);
        posP = cellfun(@(mdl) mdl.Coefficients.pValue(2), mdlPos);
        negP = cellfun(@(mdl) mdl.Coefficients.pValue(2), mdlNeg);
        tstat = cellfun(@(pos,neg) (pos.Coefficients.Estimate(2)-neg.Coefficients.Estimate(2))/norm([pos.Coefficients.SE(2) neg.Coefficients.SE(2)]), mdlPos, mdlNeg);
        diffP = arrayfun(@(t,n) tcdf(t,n), tstat, cellfun(@(D) sum(D>0),D));  % p-value of difference between the coefficients
    end
    
    % bin distances for plots
    nBins = length(edges)-1;
    D = cellfun(@(oriDiffs) sum(bsxfun(@ge,oriDiffs,edges),2), D, 'uni', false);
    
    % exclude sites that don't have enough cell pairs in each bin
    hasEnough = 200 < cellfun(@(D) min(hist(D,1:nBins)), D);
    C0 = C0(hasEnough);
    C1 = C1(hasEnough);
    S  = S(hasEnough);
    D  = D(hasEnough);
    sparsity = sparsity(hasEnough);
    fprintf('Qualifying sites n=%d\n', sum(hasEnough))
    clear hasEnough
    
    denseIx = sparsity<0.97;
    C0 = cellfun(@(C,D) accumarray(D,C,[nBins 1],@mean), C0, D, 'uni',false);
    C1 = cellfun(@(C,D) accumarray(D,C,[nBins 1],@mean), C1, D, 'uni',false);
    
    % panel B or C: noise correlations vs distance
    if vertical
        fig = Figure(1,'size',[45,32]);
        fname = fullfile(covest.plots.figPath, 'Fig6-C.eps');
    else
        fig = Figure(1,'size',[55,32]);
        fname = fullfile(covest.plots.figPath, 'Fig6-B.eps');
    end
    ticks = conv(edges,[.5 .5],'valid')';
    plot(ticks,cat(2,C0{:}),'k.-','LineWidth',0.25)
    yticks = 0:0.05:1;
    set(gca,'XTick',edges,'YTick',yticks,'YTickLabel',nozero(yticks))
    
    xlim([0 edges([end-1 end])*[.25;.75]])
    ylim([0 0.055])
    if vertical
        xlabel '\Deltaz (\mum)'
    else
        xlabel '\Deltax (\mum)'
    end
    ylabel 'avg. noise corr.'
    set(gca,'Position',[.27 .29 .67 .67])
    
    fig.cleanup
    fig.save(fname)
    
    % panel E or F: partial correlations vs distance
    if vertical
        fig = Figure(1,'size',[45,32]);
        fname = fullfile(covest.plots.figPath, 'Fig6-F.eps');
    else
        fig = Figure(1,'size',[55,32]);
        fname = fullfile(covest.plots.figPath, 'Fig6-E.eps');
    end
    ticks = conv(edges,[.5 .5],'valid')';
    plot(ticks,cat(2,C1{:}),'k.-','LineWidth',0.25);
    yticks = 0:.01:1;
    set(gca,'XTick',edges,'YTick',yticks,'YTickLabel',nozero(yticks))
    
    xlim([0 edges([end-1 end])*[.25;.75]])
    ylim([0 0.014])
    if vertical
        xlabel '\Deltaz (\mum)'
    else
        xlabel '\Deltax (\mum)'
    end
    ylabel 'avg. partial corr.'
    set(gca,'Position',[.27 .29 .67 .67])
    
    fig.cleanup
    fig.save(fname)
    
    
    clear C1 C0
    
    % panels E and F: connectivity vs distance
    S = S(denseIx);
    D = D(denseIx);
    if vertical
        fig = Figure(1,'size',[45,32]);
        fname = fullfile(covest.plots.figPath, 'Fig6-I.eps');
    else
        fig = Figure(1,'size',[55,32]);
        fname = fullfile(covest.plots.figPath, 'Fig6-H.eps');
    end
    posColor = [0 .5 0];
    negColor = [.5 0 0];
    
    neg = cellfun(@(S,D) accumarray(D,S>0,[nBins 1], @mean), S, D, 'uni',false);
    pos = cellfun(@(S,D) accumarray(D,S<0,[nBins 1], @mean), S, D, 'uni',false);
    
    plot(ticks,cat(2,pos{:}),'^-', 'Color',posColor,'MarkerSize',3,'LineWidth',.5,'MarkerFaceColor','w')
    hold on
    plot(ticks,cat(2,neg{:}), 'v:', 'Color',negColor,'MarkerSize',3,'LineWidth',.5,'MarkerFaceColor','w')
    hold off
    yticks = 0:0.1:0.2;
    set(gca,'XTick',edges,'YTick',yticks,'YTickLabel',arrayfun(@(s) strrep(sprintf('%g',s),'0.','.'), yticks,'uni',false))
    xlim([0 edges([end-1 end])*[.25;.75]])
    ylim([0 0.25])
    if vertical
        xlabel '\Deltaz (\mum)'
    else
        xlabel '\Deltax (\mum)'
    end
    ylabel 'connectivity'
    set(gca,'Position',[0.27 0.29 0.67 0.67])
    
    fig.cleanup
    fig.save(fname)
end
end

function D = distMatrix(xyz)
p = size(xyz,1);
[i,j] = ndgrid(1:p,1:p);
D = reshape(sqrt(sum((xyz(i,:)-xyz(j,:)).^2,2)),p,p);
end


function ret = offDiag(C)
% return the vector of off-diagonal elements of symmetric matrix C
p = size(C,1);
[i,j] = ndgrid(1:p,1:p);
ret = C(i<j);
end

function s=nozero(f)
s = arrayfun(@(f) strrep(sprintf('%g',f),'0.','.'), f, 'uni', false);
end