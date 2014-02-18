function demo(choice)
if nargin<1
    choice = 1;
end
close all
mag = [-0.08 0.08];  % increase size of subplots by this fraction


%% define the loss function
loss = @(S,Sigma) (trace(Sigma/S) + logDet(S))/size(S,1);   % entropy loss (aka stein's loss)
%loss = @(S,Sigma) norm(S-Sigma,'fro');                                 % sum squared error
names = {'unstructured','multifactor','sparse','sparse+lowrank'};
name =  names{choice};
letter = char(64+choice);


%% generate the true covariance matrix Sigma
rng(0)
p = 60;   % number of neurons
title_ = sprintf('truth "%s"',name);
n = ceil(p*10.0);  % number of samples in the prematrix -- regulates scale of correlations
Sigma = eye(p);
for i=1:5
    Sigma = 0.95*covest.cov(mvnrnd(zeros(p,1),Sigma,n)) + 0.05*diag(diag(Sigma));
end

switch name
    case 'unstructured'
        % do nothing
        
    case 'multifactor'
        for i=1:3
            Sigma = 0.95*covest.cov(mvnrnd(zeros(p,1),Sigma,n)) + 0.05*diag(diag(Sigma));
        end
        nlatent = 6;
        [L,Ph] = covest.factor(Sigma,nlatent);
        Sigma = L*L' + diag(Ph);
        title_ = sprintf('%s (%d factors)', title_, nlatent);
        
    case 'sparse'
        for i=1:3
            Sigma = 0.95*covest.cov(mvnrnd(zeros(p,1),Sigma,n)) + 0.05*diag(diag(Sigma));
        end
        ret = covest.lvglasso(Sigma, 0.07, 0.0, struct('refit',true,'mu',p,'max_latent',0));
        Sigma = inv(ret.S-ret.L);
        title_ = sprintf('%s (%2.1f%% sparse)', title_, 100*covest.sparsity(ret.S));
        
    case 'sparse+lowrank'
        nlatent = 3;
        for i=1:3
            Sigma = 0.95*covest.cov(mvnrnd(zeros(p,1),Sigma,n)) + 0.05*diag(diag(Sigma));
        end
        ret = covest.lvglasso(Sigma, 0.07, 0.0, struct('refit',true,'mu',p,'max_latent',nlatent));
        Sigma = inv(ret.S-ret.L);
        title_ = sprintf('%s (%2.1f%% sparse)', title_, 100*covest.sparsity(ret.S));
end

subplot 411
covest.showcov(Sigma,0.25)
title(title_)
b = get(gca, 'Position'); set(gca,'Position', b + [-b(3:4).*mag/2 b(3:4).*mag])
text(-0.3,1,char(letter+0*4),'Units','normal')

%% sample covariance Sigma0
n = 500;  % sample size
Sigma0 = covest.cov(mvnrnd(zeros(p,1),Sigma,n));
subplot 412
covest.showcov(Sigma0,0.25)
title(sprintf('sample covariance, n=%d',n))
b = get(gca, 'Position'); set(gca,'Position', b + [-b(3:4).*mag/2 b(3:4).*mag])
text(-0.3,1,char(letter+1*4),'Units','normal')


%% solve by each method, compare to truth
subplot 413
nrepeats = 10;
methods = {'sample','R.shrink','analytic shrink','shrink','unifactor','factor','glasso','lv-glasso'};
methodSelection = [1 4 6 7 8];
ns = [500 1000 2000];
colors = jet(length(methodSelection))*0.4;
styles = {'.-','.:','.--','^-.','v-'};
linewidth = [0.5 1 1 1 1 1 1 1];

if true
    for method = methodSelection
        L = nan(length(ns),nrepeats);
        fprintf('%24s',methods{method})
        for n = ns
            fprintf(' %d', n)
            for repeat = 1:nrepeats
                rng(repeat)
                X = mvnrnd(zeros(p,1),Sigma,n);
                switch methods{method}
                    case 'sample'
                        C = covest.cov(X);
                    case 'R.shrink'
                        C = covest.R_shrink(X);
                    case 'analytic shrink'
                        C = covest.shrink(X,1);
                    case 'shrink'   % cross-validated shrink
                        C = covest.crossEstimateHyper(X,loss,'shrink',1,exp(-10:0.05:0));
                    case 'unifactor'
                        C = covest.crossEstimateHyper(X,loss,'factor',1,exp(-10:0.05:0));
                    case 'factor'
                        C = covest.crossEstimateHyper(X,loss,'factor',1:20,exp(-5:0.05:0));
                    case 'glasso'
                        C = covest.crossEstimateHyper(X,loss,'lv-glasso',0,exp(-5:0.05:0));
                    case 'lv-glasso'
                        C = covest.crossEstimateHyper(X,loss,'lv-glasso',0:25,exp(-5:0.05:0));
                end
                L(ns==n,repeat) = loss(C,Sigma)-loss(Sigma,Sigma);
            end
        end
        fprintf \n
        if nrepeats >= 4
            errorbar(ns,mean(L,2),std(L,[],2),styles{method==methodSelection},...
                'LineWidth',linewidth(method==methodSelection),...
                'Color',colors(method==methodSelection,:),...
                'MarkerSize',9)
        else
            plot(ns,mean(L,2),[styles(method==methodSelection) '-'],'Color',colors(method==methodSelection,:))
        end
        hold on
        drawnow
    end
    hold off
    box off
    grid on
    legend boxoff
    legend(methods{methodSelection})
    xlim([0 ns(end)])
    xlabel('sample size','fontsize',8)
    ylabel('entropy loss','fontsize',8)
    mag = [-0.15 0.08];
    b = get(gca, 'Position'); set(gca,'Position', b + [-b(3:4).*mag/2 b(3:4).*mag])
    set(gca, 'fontsize', 6)
    text(-0.2,1,char(letter+2*4),'Units','normal')
end

% cross-validated loss
subplot 414
nrepeats = 6;
for method = methodSelection
    L = nan(length(ns),nrepeats);
    fprintf('%24s',methods{method})
    for n = ns
        fprintf(' %d', n)
        for repeat = 1:nrepeats
            rng(repeat)
            X = mvnrnd(zeros(p,1),Sigma,n);
            switch methods{method}
                case 'sample'
                    losses = covest.crossval(X,loss,'sample');
                case 'R.shrink'
                    losses = covest.crossval(X,loss,'R_shrink');
                case 'analytic shrink'
                    losses = covest.crossval(X,loss,'shrink');
                case 'shrink'   % cross-validated shrink
                    losses = covest.crossval(X,loss,'shrink',1,exp(-10:0.05:0));
                case 'unifactor'
                    losses = covest.crossval(X,loss,'factor',1,exp(-10:0.05:0));
                case 'factor'
                    losses = covest.crossval(X,loss,'factor',1:20,exp(-5:0.05:0));
                case 'glasso'
                    losses = covest.crossval(X,loss,'lv-glasso',0,exp(-5:0.05:0));
                case 'lv-glasso'
                    losses = covest.crossval(X,loss,'lv-glasso',0:25,exp(-5:0.05:0));
            end
            L(ns==n,repeat) = mean(losses);
        end
    end
    fprintf \n
    if nrepeats >= 4
        errorbar(ns,mean(L,2),std(L,[],2),styles{method==methodSelection},...
            'LineWidth',linewidth(method==methodSelection),...
            'Color',colors(method==methodSelection,:),...
            'MarkerSize',9)
    else
        plot(ns,mean(L,2),[styles(method==methodSelection) '-'],'Color',colors(method==methodSelection,:))
    end
    hold on
    drawnow
end
hold off
box off
grid on
legend boxoff
legend(methods{methodSelection})
xlim([0 ns(end)])
xlabel('sample size','fontsize',8)
ylabel('entropy loss','fontsize',8)
mag = [-0.15 0.08];
b = get(gca, 'Position'); set(gca,'Position', b + [-b(3:4).*mag/2 b(3:4).*mag])
set(gca, 'fontsize', 6)
text(-0.2,1,char(letter+2*4),'Units','normal')



set(gcf,'PaperSize',[3 8], 'PaperPosition', [0 0 3 8])
print('-dtiff','-r900',sprintf('fig2-%s',letter))