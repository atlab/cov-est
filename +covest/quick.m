function quick

c = covest.CovMatrix & 'nfolds>1';

c1 = pro(c, 'method->m1','cv_loss->l1');
c2 = pro(c, 'method->m2','cv_loss->l2');

r = 'm1=70 and m2=90';
%[l1,l2] = fetchn(c1*c2 & r, 'l1', 'l2');
[l1,l2] = fetchn(covest.ActiveCells, c1*c2 & r, 'avg(l1)->ll1','avg(l2)->ll2');
l1(isnan(l1))=inf;
fprintf('Median %g, p=%1.0e, frac positive %g\n', ...
    median(l1-l2), signrank(l1,l2), mean(l1>l2))
hist(l1-l2,30)
figure(gcf)
end