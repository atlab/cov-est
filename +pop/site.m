function site(scanNum,isDepth)

if isDepth
    bounds = [0 30 60];
    dZ = 1000;
    dX = 50;
else
    bounds = [0 25 75 125];
    dZ = 50;
    dX = 1000;
end
key = fetch(pop.AodBinnedTraces & sprintf('mod(aod_scan_start_time,1e4)=%d',scanNum) & 'bin_opt=0');

% covariance estimate
cove = fetch(...
    pop.AodCovEstimate*pro(pop.AodBinnedTraces,'cellnums') ...
    & key & 'cov_estim_num=14', '*');
[x,y,z,xyzcellnums] = fetchn(aod.Traces*aod.UniqueCell & key,'x','y','z','cell_num');
Sigma0 = fetch1(pop.AodCovEstimate & key & 'cov_estim_num=0','cov_matrix');
C0 = corrcov(Sigma0);
p = size(C0,1);

% deltaX, deltaZ matrices
deltaX = nan(size(cove.sparse));
deltaZ = nan(size(cove.sparse));
nCells = size(cove.sparse,1);
for i=2:nCells
    ic = find(xyzcellnums==cove.cellnums(i));
    assert(length(ic)==1)
    for j=1:i-1
        jc = find(xyzcellnums==cove.cellnums(j));
        assert(length(jc)==1)
        
        d = norm([x(ic)-x(jc) y(ic)-y(jc)],2);
        deltaX(i,j) = d;
        deltaX(j,i) = d;
        
        d = abs(z(ic)-z(jc));
        deltaZ(i,j) = d;
        deltaZ(j,i) = d;
    end
end
sz = size(cove.sparse);
[ix,jx] = ndgrid(1:sz(1),1:sz(2));
restrict = deltaZ < dZ & deltaX < dX;              
offDiag = jx>ix & restrict;
offDiag = offDiag(:);
connected = offDiag & cove.sparse(:);
thresholded = offDiag & abs(C0(:)) > quantile(abs(C0(offDiag)),cove.sparsity);
rec.total.deltaZ = deltaZ(offDiag);
rec.total.deltaX = deltaX(offDiag);
rec.connected.deltaZ = deltaZ(connected);
rec.connected.deltaX = deltaX(connected);
rec.thresh.deltaZ = deltaZ(thresholded);
rec.thresh.deltaX = deltaX(thresholded);

if isDepth
    use.total = rec.total.deltaZ;
    use.connected = rec.connected.deltaZ;
    use.thresh = rec.thresh.deltaZ;
else
    use.total = rec.total.deltaX;
    use.connected = rec.connected.deltaX;
    use.thresh = rec.thresh.deltaX;
end
clear rec

use.C0 = C0(offDiag);

fprintf('%d cells, sparsity %2.1f%%, mean connectivity %2.2f per cell\n', p, 100*cove.sparsity, sum(~~cove.sparse(jx(:)~=ix(:)))/p)
fprintf('overlap %2.1f%%\n',100*sum(thresholded.*connected)/sum(connected))


bounds = [bounds max(use.total)];
meanAbsCorrs = arrayfun(@(i) mean(abs(use.C0(between(use.total,bounds(i+(-1:0)))))), 2:length(bounds));
edgesPerBin = arrayfun(@(i) sum(between(use.connected,bounds(i+(-1:0)))), 2:length(bounds))
connProb = arrayfun(@(i) ...
    sum(between(use.connected,bounds(i+(-1:0))))/...
    sum(between(use.total,bounds(i+(-1:0)))), 2:length(bounds));
threshProb = arrayfun(@(i) ...
    sum(between(use.thresh,bounds(i+(-1:0))))/...
    sum(between(use.total,bounds(i+(-1:0)))), 2:length(bounds));
x= 0.5:length(bounds)-1;
c1 = [0.5 0.3 0];
c2 = [0 0.5 0.3];
[ax,h1,h2] = plotyy(x,meanAbsCorrs,x,connProb);
hold(ax(2),'on'),
h3 = plot(ax(2),x,threshProb,'LineStyle',':','Marker','v','LineWidth',2,'MarkerSize',9);
hold(ax(2),'off')
set(h1,'Color',c1,'LineWidth',1,'LineStyle','-','Marker','s','MarkerSize',9)
set(h2,'Color',c2,'LineWidth',2,'LineStyle','-','Marker','o','MarkerSize',9)
set(ax,'XTick',0:length(bounds)-2,'XTickLabel',bounds(1:end-1),'YTick',0:.01:0.1,'XLim',[0 length(bounds)-0.5],'box','off','XGrid','on','YGrid','on','Ylim',[0 0.08]/scale)
set(ax(1),'YColor',c1)
set(ax(2),'YColor',c2)
set(ax,'Position',get(ax(1),'Position').*[1 1 0.95 1])
if isDepth
    xlabel 'vertical distance (\mum)'
else
    xlabel 'lateral distance (\mum)'
end
ylabel(ax(1),'avg abs correlation')
ylabel(ax(2),'connectivity')
pos = get(ax(1),'Position');
set(ax,'Position',pos.*[1 1.2 1 1])
%legend('avg corr','estimator D','thresh corr')
%h = legend('location','NorthEastOutside')
%set(h,'fontsize',15)


%set(gcf,'PaperSize',[4 3],'PaperPosition',[0 0 4 3])
%print -dpdf thresh_vs_distance
end


function yes = between(a,rng)
yes = a > rng(1) & a<=rng(2);
end

