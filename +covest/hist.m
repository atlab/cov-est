function hist(C)
% histogram of correlations
[x,y] = meshgrid(1:size(C,2),1:size(C,1));
C = corrcov(C);
C = C(x<y);
m = max(abs(C));
ix = linspace(-m,m,ceil(sqrt(numel(C)/2)));
hist(C,ix)
grid on
box off
xlabel correlation
hold on
plot([0 0]+mean(C), ylim, 'r')
plot([0 0], ylim, 'k')
plot([0 0]+mean(C)+std(C), ylim, 'b')
plot([0 0]+mean(C)-std(C), ylim, 'b')
hold off
colormap((1-jet)/4+0.2)
xlim(ix([1 end]))