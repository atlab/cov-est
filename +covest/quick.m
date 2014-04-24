function quick

c = covest.CovMatrix & 'nfolds>1';
c = c.pro('cv_loss->cl');
%c = pro(covest.Method*covest.ActiveCells, c, 'avg(cv_loss)->cl');

c1 = pro(c & 'method=10', 'method->m1','cl->l1');
c2 = pro(c & 'method=50', 'method->m2','cl->l2');

[l1,l2] = fetchn(c1*c2, 'l1', 'l2');
l1(isnan(l1))=inf;
median(l1-l2)
signrank(l1-l2)
hist(l1-l2,100)
figure(gcf)
end