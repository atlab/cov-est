function supp1A
letter = 'A';
for alpha = [.01 .05 .1 .5 .6 .65 .7 .75 .8 .85 .9]
    nfolds = 10;
    L1 = [];
    L2 = [];
    for key = fetch(covest.ActiveCells & 'preprocess_method_num=5' & 'high_repeats' & 'ncells>100')'
        [X,evokedBins,sel] = fetch1(covest.ActiveCells*covest.Traces & key, ...
            'trace_segments','evoked_bins', 'selection');
        X = X(1:min(end,evokedBins),:,:,sel);
        [nBins,nConds,nTrials,nCells] = size(X);
        
        L1(end+1) = mean(arrayfun(@(k) getLoss(X,k,true), 1:nfolds)); %#ok<AGROW>
        L2(end+1) = mean(arrayfun(@(k) getLoss(X,k,false), 1:nfolds));     %#ok<AGROW>
    end
    fprintf('Median difference = %3.5g, p-value=%1.1e\n', median(L1-L2), signrank(L1-L2))
fig = Figure(1, 'size', [83 20]);
h = boxplot(L1-L2,'jitter',0,'colors','k',...
    'labels',sprintf('a = %1.2f',alpha),'orientation','horizontal','outliersize',3);
set(h(1:2,:),'LineStyle','-','LineWidth',.25)
set(h(7,:),'MarkerEdgeColor','k')
xlabel nats/cell/bin
set(gca,'YColor',[1 1 1],'YDir','reverse')
hold on
plot([0 0],ylim,'k:')
hold off
axis tight
set(gca,'Position',[.20 .5 0.75 0.5])
fig.cleanup
fig.save(fullfile(covest.plots.figPath, sprintf('Supp1-%c-bin.eps',letter)))
letter = letter+1;
end



    function L = getLoss(X, k, binwise)
        loss = @(S,Sigma)(trace(Sigma/S)+cove.logDet(S))/size(S,1);

        [Z, Y] = cove.splitTrials(X,k,nfolds);
        Z = double(Z);
        
        % training
        M = nanmean(Z,3);
        Z = bsxfun(@minus, Z, M);
        % common variance
        V0 = nanvar(reshape(Z,[],nCells));
        V0 = reshape(V0,1,1,1,nCells);
        if ~binwise
            V = V0;
        else
            % bin-wise variance
            V = nanvar(Z,[],3);
            % condition-wise variance
            % V = mean(V);
            % regularize: bias toward common variance
            V = bsxfun(@plus,(1-alpha)*V, alpha*V0);
        end
        V = sqrt(V);
        Z = bsxfun(@rdivide, Z, V);  % zscore
        C = corrcov(cove.cov(reshape(Z,[],nCells)));  % correlation matrix
        
        % testing
        Y = double(Y);
        Y = bsxfun(@minus, Y, M);
        Y = bsxfun(@rdivide, Y, V);
        
        L = loss(C, cove.cov(reshape(Y,[],nCells)));        
    end
end