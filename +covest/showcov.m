function showcov(C,range,plotBar)
% display the covariance matrix C
if nargin<2
    range = (max(max(abs(C))));
end
if nargin<3
    plotBar = true;
end
imagesc(C,[-1 1]*range)
axis image
colormap(covest.doppler)
if plotBar
    colorbar
end
axis off
set(gca,'XTick',[],'YTick',[])
set(gca,'fontsize',8)