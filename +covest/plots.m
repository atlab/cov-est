% plots and figures for the covariance estimation paper


classdef plots
    
    properties(Constant)
        figPath = '~/cov/figures/src/'
        exampleSite = 'mod(aod_scan_start_time,10000)=4328 && high_repeats'
    end
    
    methods(Static)
        
        
        function supp2pre
            loss = @(S,Sigma)(trace(Sigma/S)+cove.logDet(S))/size(S,1);
            
            
            % X := nBins * nDirs * nTrials * nCells
            [X,selection] = fetch1(covest.Traces*covest.ActiveCells & covest.plots.exampleSite, ...
                'trace_segments', 'selection');
            X = double(X(:,:,:,selection));
            Xsaved = X;
            
            [hypers, evokedBins] = fetch1(...
                covest.Traces*covest.CovMatrix & covest.plots.exampleSite ...
                & 'method=90' & 'nFolds = 1','hypers','evoked_bins');
            
            nFolds = 10;
            cvloss = nan(41,41,nFolds);
            [alpha,beta] = ndgrid(...
                hypers(2)*exp(linspace(-1,1,size(cvloss,1))), ...
                hypers(3)*exp(linspace(-1,1,size(cvloss,2))));
            sparsity = nan(size(alpha));
            nLatent = nan(size(alpha));
            for iAlpha = 1:size(cvloss,1)
                for iBeta = 1:size(cvloss,2)
                    hypers(2) = alpha(iAlpha,iBeta);
                    hypers(3) = beta(iAlpha,iBeta);
                    fprintf('%02d-%02d ',iAlpha,iBeta)
                    for k=1:nFolds
                        [X, XTest] = cove.splitTrials(Xsaved,k,nFolds);
                        [C,M,V] = cove.estimate(X, [], [], evokedBins, 'lv-glasso', hypers);
                        CTest = cove.estimate(XTest,M,V,evokedBins, 'sample', []);
                        cvloss(iAlpha,iBeta,k) = loss(C,CTest);
                        fprintf .
                    end
                    [~,~,extras] = cove.estimate(Xsaved, [], [], evokedBins, 'lv-glasso', hypers);
                    sparsity(iAlpha,iBeta) = cove.sparsity(extras.S);
                    nLatent(iAlpha,iBeta) = size(extras.H,2);
                    fprintf(' sparsity %1.3f  nLatent %3d\n', ...
                        sparsity(iAlpha,iBeta), nLatent(iAlpha,iBeta));
                end
            end
            save ~/comment4v4
        end
        
        
        function supp2
            s = load('~/comment4v5');
            connectivity = 1-s.sparsity;
            
            [hypers] = fetch1(...
                covest.Traces*covest.CovMatrix & covest.plots.exampleSite ...
                & 'method=90' & 'nFolds = 1','hypers');
            fig = Figure(1,'size',[80 70]);
            cvloss = mean(s.cvloss,3);
            cvloss = cvloss - min(cvloss(:));
            alpha = s.alpha(:,1);
            beta = s.beta(1,:);
            [C,h] = contour(alpha,beta,cvloss','k');
            clabel(C,h)
            xlabel \alpha
            ylabel \beta
            yticks = 0.1:0.1:1;
            xticks = 0.01:0.01:0.1;
            set(gca,'XScale','log','YScale','log',...
                'XTick',xticks,'XTickLabel',nozero(xticks), ...
                'YTick',yticks,'YTickLabel',nozero(yticks))
            hold on, plot(hypers(2),hypers(3),'r+','MarkerSize',30), hold off
            fig.cleanup
            fig.save(fullfile(covest.plots.figPath,'Supp2-A.eps'))
            
            fig = Figure(1,'size',[80 70]);
            
            v = 0:0.01:0.5;
            vv = v;
            vv(1:5:end) = [];
            v = v(1:5:end);
            contour(alpha,beta,connectivity',vv,'Color',[1 1 1]*.7,'LineWidth',.25)
            hold on
            [C,h] = contour(alpha,beta,connectivity',v,'k');
            hold off
            clabel(C,h)
            xlabel \alpha
            ylabel \beta
            yticks = 0.1:0.1:1;
            xticks = 0.01:0.01:0.1;
            set(gca,'XScale','log','YScale','log',...
                'XTick',xticks,'XTickLabel',nozero(xticks), ...
                'YTick',yticks,'YTickLabel',nozero(yticks))
            hold on, plot(hypers(2),hypers(3),'r+','MarkerSize',30), hold off
            fig.cleanup
            fig.save(fullfile(covest.plots.figPath,'Supp2-B.eps'))
            
            
            fig = Figure(1,'size',[80 70]);
            [C,h] = contour(alpha,beta,s.nLatent','k');
            clabel(C,h)
            xlabel \alpha
            ylabel \beta
            yticks = 0.1:0.1:1;
            xticks = 0.01:0.01:0.1;
            set(gca,'XScale','log','YScale','log',...
                'XTick',xticks,'XTickLabel',nozero(xticks), ...
                'YTick',yticks,'YTickLabel',nozero(yticks))
            hold on, plot(hypers(2),hypers(3),'r+','MarkerSize',30), hold off
            fig.cleanup
            fig.save(fullfile(covest.plots.figPath,'Supp2-C.eps'))
            
            
            fig = Figure(1,'size',[80 70]);
            [n,sp] = fetch1(covest.CovMatrix & covest.plots.exampleSite ...
                & 'method=90' & 'nFolds = 1','lowrank','sparsity');
            n = size(n,2);
            F = scatteredInterpolant(connectivity(:),s.nLatent(:),cvloss(:),'linear','none');
            [nLatent,connectivity] = ndgrid(0:1:120,0.0:0.01:0.3);
            cvloss = F(connectivity,nLatent);
            cvloss = medfilt1(cvloss,7);
            k = hamming(3); k = k/sum(k);
            cvloss = imfilter(imfilter(cvloss,k,nan),k',nan);
            [C,h] = contour(connectivity,nLatent,cvloss,'k');
            clabel(C,h)
            hold on, plot(1-sp,n,'r+','MarkerSize',30), hold off
            xlabel connectivity
            ylabel '# latent'
            ylim([0 110])
            fig.cleanup
            fig.save(fullfile(covest.plots.figPath,'Supp2-D.eps'))
        end
        
        
        function stimCond
            % correlations and connectivity conditioned on the stimulus
            doCorrs = false;
            
            c = covest.CovMatrix & 'nfolds=1';
            c1 = c.pro('method->m1','sparse->s1','cov_matrix->c1');
            c2 = c.pro('method->m2','sparse->s2','cov_matrix->c2');
            [S1,S2,C1,C2,dirs,pval,selection] = fetchn(covest.ActiveCells*pro(covest.OriTuning,'high_repeats->temp','von_pref','von_p_value')*c1*c2 & 'm1=92 && m2=93', ...
                's2','s1','c1','c2','von_pref','von_p_value','selection');
            dirs = cellfun(@(dirs,selection) dirs(selection)/pi*180, dirs, selection, 'uni', false);
            pval = cellfun(@(pval,selection) pval(selection), pval, selection, 'uni', false);
            if ~doCorrs
                d =[
                    cellfun(@(S,dirs,pval) ndegree(S,dirs,pval,90), S1, dirs, pval)'
                    cellfun(@(S,dirs,pval) ndegree(S,dirs,pval, 0), S1, dirs, pval)'
                    cellfun(@(S,dirs,pval) ndegree(S,dirs,pval,90), S2, dirs, pval)'
                    cellfun(@(S,dirs,pval) ndegree(S,dirs,pval, 0), S2, dirs, pval)'
                    ];
            else
                d =[
                    cellfun(@(C,dirs,pval) avgCorr(C,dirs,pval,90), C1, dirs, pval)'
                    cellfun(@(C,dirs,pval) avgCorr(C,dirs,pval, 0), C1, dirs, pval)'
                    cellfun(@(C,dirs,pval) avgCorr(C,dirs,pval,90), C2, dirs, pval)'
                    cellfun(@(C,dirs,pval) avgCorr(C,dirs,pval, 0), C2, dirs, pval)'
                    ];
            end
            
            d = d(:,~any(isnan(d)));
            fig = Figure(1,'size',[40 40]);
            plot(d([1,3],:),'-ro')
            hold on
            plot(d([2,4],:),'-gx')
            xlim([0.5 2.5])
            fig.cleanup
            if doCorrs
                fig.save(fullfile(covest.plots.figPath,'Supp5-A.eps'))
            else
                fig.save(fullfile(covest.plots.figPath,'Supp5-B.eps'))
            end
            
            function ret = ndegree(S,dirs,pval,filt)
                ix = pval < 0.01 & oriDiff(filt,dirs') < 40;
                if sum(ix)<10
                    ret = nan;
                else
                    ret = (1-cove.sparsity(S(ix,ix)))/(1-cove.sparsity(S));
                end
            end
            
            function ret = avgCorr(C,dirs,pval,filt)
                ix = pval < 0.01 & oriDiff(filt,dirs') < 40;
                if sum(ix)<10
                    ret = nan;
                else
                    ret = cove.avgCorr(C(ix,ix));%/cove.avgCorr(C);
                end
            end
            
        end
        
        
        function fig2
            % panel C
            fig = Figure(1,'size',[35 32]);
            
            C = fetch1(covest.CovMatrix & covest.plots.exampleSite & 'method=0' & 'nfolds=1','cov_matrix');
            p = size(C,1);
            C = corrcov(C);
            imagesc(C,.1*[-1 1]);
            colormap(cove.doppler)
            axis image
            set(gca,'YTick',[1 p],'XTick',[])
            set(gca, 'Position', [0.20 0.04 0.795 0.89]);
            fig.cleanup
            box on
            fig.save(fullfile(covest.plots.figPath,'Fig2-E.eps'))
            
            fig = Figure(1,'size',[35 30]);
            [x,y] = meshgrid(1:p,1:p);
            [a,b] = hist(C(x<y),-0.6:0.002:0.2);
            area(b,a,'LineWidth',0.5,'EdgeColor',[0 0 0.4],'FaceColor',[0.7 0.7 0.8]);
            hold on
            plot([0 0],[0 1]*max(a)*1.1,'k')
            plot([1 1]*mean(C(x<y)),[0 1]*max(a)*1.1,'r-','LineWidth',1)
            xlim([-0.06 0.21])
            ylim([0 1.1*max(a)])
            xlabel correlations
            set(gca,'YTick',[])
            ticks = -.05:0.05:0.2;
            set(gca,'XTick',ticks,'XTickLabel',nozero(ticks))
            set(gca,'Position',[0.1 .3 1.0 0.7])
            set(gca,'YColor',[1 1 1]*.999)
            fig.cleanup
            fig.save(fullfile(covest.plots.figPath,'Fig2-F.eps'))
        end
        
        
        
        function fig3
            for f = {'Fig3','Supp3','Supp4'}
                switch f{1}
                    case {'Fig3','Supp4'}
                        pairs = {
                            0   90    'sample'
                            10  90    'diag'
                            30  90    'factor'
                            80  90    'sparse'
                            };
                    case 'Supp3'
                        pairs = {
                            0   80    'sample'
                            10  80    'diag'
                            30  80    'factor'
                            90  80    sprintf('sparse\n+latent')
                            };
                    otherwise
                        error 'unknown figure'
                end
                ticks = 0:0.01:.1;
                if any(strcmp(f{1},{'Fig3','Supp3'}))
                    c = covest.CovMatrix & 'nfolds>1';
                    c1 = pro(c, 'method->m1','cv_loss->l1');
                    c2 = pro(c, 'method->m2','cv_loss->l2');
                    xlabl = 'nats/cell/bin';
                elseif strcmp(f{1},'Supp4')
                    c = covest.CovMatrix*covest.QuadLoss & 'nfolds>1';
                    c1 = pro(c, 'method->m1','quad_loss->l1');
                    c2 = pro(c, 'method->m2','quad_loss->l2');
                    ticks = (0:0.01:.1)/20;
                    xlabl = 'quadratic loss';
                else
                    error 'unknown figure'
                end
                x = arrayfun(@(i) ...
                    fetchn(covest.ActiveCells, c1*c2 & sprintf('m1=%d and m2=%d',pairs{i,1:2}),...
                    'avg(l1)-avg(l2)->diff'),...
                    1:size(pairs,1), 'uni', false);
                x = [x{:}];
                
                fprintf('Medians:%s\n', sprintf(' %1.2e',median(x)))
                p = arrayfun(@(i) signrank(x(:,i)), 1:size(x,2));
                fprintf('p-values: %s\n', sprintf(' %1.1e',p))
                fig = Figure(1, 'size', [163 35]);
                h = boxplot(x,'jitter',0,'colors','k',...
                    'labels',pairs(:,3),'orientation','horizontal','outliersize',3);
                set(h(1:2,:),'LineStyle','-','LineWidth',.25)
                set(h(7,:),'MarkerEdgeColor','k')
                xlabel(xlabl)
                set(gca,'XTick',ticks,'XTickLabel',nozero(ticks))
                set(gca,'YColor',[1 1 1],'YDir','reverse')
                hold on
                plot([0 0],ylim,'k:')
                hold off
                axis tight
                set(gca,'Position',[.08 .3 0.90 0.7])
                fig.cleanup
                fig.save(fullfile(covest.plots.figPath, [f{1} '.eps']))
            end
        end
        
        
        function fig4
            select = [covest.plots.exampleSite ' && nfolds=1'];
            C0 = fetch1(covest.CovMatrix & select & 'method=0','corr_matrix');
            [C1,S,L] = fetch1(covest.CovMatrix & select & 'method=90','corr_matrix','sparse','lowrank');
            p = size(C0,1);
            CC0 = corrcov(C0);
            
            fprintf('sparsity: %2.1f%%, avg node degree = %3.2f, low-rank=%d\n', ...
                cove.sparsity(S)*100, cove.nodeDegree(S), size(L,2))
            
            % correlation matrix: sample / regularized
            comboPlot(CC0,corrcov(C1),fullfile(covest.plots.figPath, 'Fig4-A.eps'))
            densityPlot(corrcov(C0),corrcov(C1),fullfile(covest.plots.figPath,'Fig4-B'))
            
            % partial correlation matrix: sample / regularized
            iC0 = -corrcov(inv(C0));
            iC0 = (~eye(p)).*iC0;
            iC = inv(C1);
            ds = diag(sqrt(diag(iC)));
            %ds = diag(sqrt(diag(S)));
            iC = -(~eye(p)).*(ds\iC/ds);
            
            % partial correlation matrix: sample / regularized
            comboPlot(iC0,iC,fullfile(covest.plots.figPath, 'Fig4-C.eps'))
            densityPlot(iC0,iC,fullfile(covest.plots.figPath, 'Fig4-D'))
            % replace interaction matrix with thresholded correlations
            
            comboPlot(-ds\S/ds,ds\L*L'/ds,...
                fullfile(covest.plots.figPath, 'Fig4-E.eps'))
            densityPlot(CC0,-ds\S/ds,...
                fullfile(covest.plots.figPath, 'Fig4-F'),[0 0.2],true)
            
            
            
            function densityPlot(C1,C2,filename,axisAlphas,threshold)
                if nargin<4
                    axisAlphas=[.2 .2];
                end
                rng = 0.1;
                ctrs = linspace(-rng*0.8,rng*1.8,200);
                pp = size(C1,1);
                [i,j] = ndgrid(1:pp,1:pp);
                I = hist3([C2(i<j) C1(i<j)],{ctrs ctrs});
                cmap=[[1 1 1];jet(100)];
                I = reshape(cmap(min(end,I+1),:), [size(I) 3]);
                if nargin>=5
                    threshold = quantile(abs(C1(i>j)),cove.sparsity(C2));
                    temp = I(:,abs(ctrs)<threshold,:);
                    temp(repmat(all(temp==1,3),1,1,3)) = .8;
                    I(:,abs(ctrs)<threshold,:) = temp;
                end
                image(ctrs,ctrs,I)
                set(gca,'XTick',rng,'YTick',rng)
                set(refline(1,0),'Color','r','LineWidth',.5)
                set(gca,'Position',[0.04 0.04 0.92 0.92])
                axis xy
                axis square
                hold on
                PlotAxisAtOrigin(axisAlphas)
                hold off
                set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 3.5 3.5],...
                    'PaperSize',[3.5 3.5])
                print('-dpng','-r800',filename)
            end
            
            function comboPlot(C1,C2,filename)
                fig = Figure(1,'size',[50 46]);
                d = 8;  % zoom factor for plotting matrices
                pp = size(C1,1);
                [i,j] = ndgrid(1:pp,1:pp);
                gap = 15;
                C = zeros(size(C1)+[0 gap]);
                C(:,1:pp) = C(:,1:pp) + C1.*(i>j);
                C(:,(1:pp)+gap) = C(:,(1:pp)+gap) + C2.*(i<j);
                imagesc(0:pp+gap+1,0:pp+1,...
                    imresize(C,d,'nearest'), .1*[-1 1])
                axis image
                colormap(cove.doppler)
                set(gca,'YTick',[1 pp],'XTick',[])
                pos = [.14 .03 .85 .945];
                set(gca,'Position',pos)
                line([0 pp+1]+.5,[0 pp+1],'Color','k','LineWidth',.25)
                line([0 pp+1]+gap-.5,[0 pp+1],'Color','k','LineWidth',.25)
                fig.cleanup
                box on
                fig.save(filename)
            end
            
        end
        
        
        function network
            % figures 4E, 4G, H,a and I
            
            clf
            
            doInteractions = true;  % true=plot interactions, false=only cells
            doCorr = true; % when true, plot thresholded correlations in the fragment
            doFragment = true;
            doFragment = doFragment || doCorr;
            
            % figure 4-G,H,I
            alpha = 0.05;  % tuning signficance levels
            zref = 200;  % cortical depth of the center of the scan
            xoffset = -30;
            yoffset = -26;
            zoffset = 22;
            
            zticks = 100:50:300;
            xticks = -200:50:200;
            yticks = -200:50:200;
            zm = .60;
            panx = 1.3;
            pany = 0;
            alphaMultiplier = 3.0+3.0*doFragment;
            
            fname = fullfile(covest.plots.figPath,'Fig4-G');
            
            if doFragment
                fragmentRadius = 48;
                paperSize = [6.5 5.0];
                lineWidth = 1;
                zm = 0.7;
                panx = -2.5;
            else
                fragmentRadius = inf;
                if doInteractions
                    paperSize = [12 11];
                    lineWidth = 0.5;
                else
                    panx = 2.5;
                    paperSize = [8.3 7.5];
                end
            end
            
            % get cell positions, tuning, and sparse interactions
            key = fetch(covest.CovMatrix & covest.plots.exampleSite & 'method=90' & 'nfolds=1');
            assert(isscalar(key))
            xyz = fetch1(covest.Traces & key & 'high_repeats','cell_xyz');
            selection = fetch1(covest.ActiveCells & key, 'selection');
            [ori,pval] = fetch1(covest.OriTuning & rmfield(key,'high_repeats'), 'von_pref','von_p_value');
            [S,L] = fetch1(covest.CovMatrix & key, 'sparse', 'lowrank');
            S = -corrcov(S);  % convert to partial correlations
            
            % thresholded correlations
            C0= corrcov(fetch1(covest.CovMatrix & setfield(key, 'method', 0), 'corr_matrix')); %#ok<SFLD>
            sparsity = fetch1(covest.CovMatrix & key, 'sparsity');
            p = size(C0,1);
            [i,j] = ndgrid(1:p,1:p);
            C0 = C0.*(abs(C0)>quantile(abs(C0(i<j)), sparsity));
            
            % report everything
            fprintf('Ssparsity = %2.1f%%\n', 100*sparsity)
            fprintf('Overlap = %2.1f%%\n',  100*sum(C0(i<j) & S(i<j))/sum(~~S(i<j)))
            fprintf('Latent = %d\n', size(L,2));
            fprintf('Negative interactions = %2.1f%%\n', 100*sum(S(i<j)<0)/sum(~~S(i<j)))
            
            if doInteractions
                if doCorr
                    S = C0;
                    fname = fullfile(covest.plots.figPath,'Fig4-I');
                elseif doFragment
                    fname = fullfile(covest.plots.figPath,'Fig4-H');
                end
            else
                fname = fullfile(covest.plots.figPath,'Fig2-E');
            end
            
            x = xyz(:,1);
            y = xyz(:,2);
            z = xyz(:,3)+zref;
            
            % plot balls
            hue = mod(ori(:)/pi,1);
            sat = pval(:)<alpha;
            val = (1-.2*(pval(:)>=alpha)).*selection(:);
            color = hsv2rgb([hue sat val]);
            clear hue sat val
            
            fragIdx = (x-xoffset).^2+(y-yoffset).^2+(z-zref-zoffset).^2 < fragmentRadius^2;
            
            scatter3sph(x(fragIdx),y(fragIdx),z(fragIdx),'siz',3,'col',color(fragIdx,:))
            light('Position',[0.5 0.5 1],'Style','infinit','Color',[1 1 1])
            lighting gouraud
            axis vis3d
            axis equal
            set(gca,'ZDir','reverse')
            camproj perspective
            grid on
            
            % show interactions
            x = x(selection);
            y = y(selection);
            z = z(selection);
            
            if doFragment
                fragIdx = (x-xoffset).^2+(y-yoffset).^2+(z-zref-zoffset).^2 < fragmentRadius^2;
                x = x(fragIdx);
                y = y(fragIdx);
                z = z(fragIdx);
                S = S(fragIdx,fragIdx);
            end
            
            
            if doInteractions
                [i,j] = ndgrid(1:size(S,1),1:size(S,2));
                % positive interactions
                ix = find(j(:)>i(:) & S(:)>0);
                for ix = ix'
                    lineAlpha = min(1,(alphaMultiplier*2.0*abs(S(ix))));
                    patch([x(i(ix)) x(j(ix))],[y(i(ix)) y(j(ix))], [z(i(ix)) z(j(ix))], [0 0], ...
                        'EdgeColor',[0 1 0]/2,'EdgeAlpha',lineAlpha,'LineWidth',lineWidth,'FaceColor','none')
                end
                
                % negative interactions
                ix = find(j(:)>i(:) & S(:)<0);
                for ix = ix'
                    lineAlpha = min(1,(alphaMultiplier*4.0*abs(S(ix))));
                    patch([x(i(ix)) x(j(ix))],[y(i(ix)) y(j(ix))], [z(i(ix)) z(j(ix))], [0 0], ...
                        'EdgeColor',[1 0 0]/2,'EdgeAlpha',lineAlpha,'LineWidth',lineWidth,'FaceColor','none')
                end
            end
            view(25-90, 65)
            zoom(zm)
            camPos = get(gca,'CameraPosition');
            set(gca,'CameraPosition',camPos*0.5)
            set(gca,'ZTick',zticks)
            set(gca,'XTick',xticks,'YTick',yticks)
            campan(panx,pany)
            set(gca,'fontsize',8,'linewidth',0.25,'TickLength',get(gca,'TickLength')*0.75)
            
            set(gcf,'PaperUnits','centimeters','PaperSize',paperSize,'PaperPosition',[0 0 paperSize])
            print('-dpdf','-r800',fname)
        end
        
        
        function fig5
            % panel A:  #latent vs #neurons
            [L,S,nCells,highlight] = fetchn(covest.CovMatrix*covest.ActiveCells & 'method=90' & 'nfolds=1', ...
                'lowrank', 'sparse', 'ncells', '(mod(aod_scan_start_time,10000)=4328)->highlight');
            fname = fullfile(covest.plots.figPath, 'Fig5-A.eps');
            fig = Figure(1,'size',[40 40]);
            
            nLatent = cellfun(@(L) size(L,2), L);
            scatter(nCells,nLatent,17,nCells,'filled')
            colormap(flipud(hot)/1.5)
            hix = find(highlight);
            hold on, plot(nCells(hix),nLatent(hix),'bo')
            set(gca,'XTick',[0 150 300])
            xlabel '# neurons'
            set(gca,'YTick',0:25:100)
            ylabel '# latent'
            
            fig.cleanup
            fig.save(fname)
            
            % panel B:  #node degree vs #neurons
            fname = fullfile(covest.plots.figPath, 'Fig5-B.eps');
            fig = Figure(1,'size',[40 40]);
            
            nodeDegree = cellfun(@cove.nodeDegree, S);
            scatter(nCells,nodeDegree,17,nCells,'filled')
            colormap(flipud(hot)/1.5)
            hold on, plot(nCells(hix),nodeDegree(hix),'bo')
            set(gca,'XTick',[0 150 300])
            xlabel '# neurons'
            set(gca,'YTick',0:25:100)
            ylabel '# avg node degree'
            
            fig.cleanup
            fig.save(fname)
            
            
            % Panel C: avg node degree vs # latent
            fname = fullfile(covest.plots.figPath, 'Fig5-C.eps');
            fig = Figure(1,'size',[40 40]);
            
            scatter(nLatent,nodeDegree,17,nCells,'filled')
            colormap(flipud(hot)/1.5)
            hold on, plot(nLatent(hix),nodeDegree(hix),'bo')
            set(gca,'XTick',0:25:100)
            xlabel '# latent'
            set(gca,'YTick',0:25:100)
            ylabel '# avg node degree'
            
            fig.cleanup
            fig.save(fname)
            
            
            % Panel D: average regularized partial correlation vs avg sample correlation
            fname = fullfile(covest.plots.figPath, 'Fig5-D.eps');
            fig = Figure(1,'size',[40 40]);
            c = covest.CovMatrix & 'nfolds=1';
            c0 = c.pro('method->m0','cov_matrix->c1');
            c1 = c.pro('method->m1','cov_matrix->c2','sparse');
            [C0,C1,S1,hix2] = fetchn(c0*c1 & 'm0=0' & 'm1=90', ...
                'c1', 'c2', 'sparse', '(mod(aod_scan_start_time,10000)=4328)->highlight');
            hix2 = find(hix2);
            m0 = cellfun(@cove.avgCorr, C0);
            m1 = -cellfun(@(C) cove.avgCorr(inv(C)), C1);
            
            scatter(m0,m1,17,nCells,'filled')
            colormap(flipud(hot)/1.5)
            hold on, plot(m0(hix2),m1(hix2),'bo')
            ticks = 0:.02:.1; set(gca,'XTick',ticks,'XTickLabel',nozero(ticks))
            xlabel 'avg sample corr'
            ticks = 0:.002:.01; set(gca,'YTick',ticks,'YTickLabel',nozero(ticks))
            ylabel '# avg reg partial corr'
            ylim(ylim.*[0 1])
            xlim(xlim.*[0 1])
            fprintf('Averages: corrs %1.2e, pcorrs %1.2e\n', ...
                mean(m0), mean(m1))
            fprintf('Coefficients of variation:  corrs %1.2f, pcorrs %1.2f\n', ...
                std(m0)/mean(m0), std(m1)/mean(m1))
            
            fig.cleanup
            fig.save(fname)
            
            
            % Panel E:  % overlap with high correlations vs connectivity
            fname = fullfile(covest.plots.figPath, 'Fig5-E.eps');
            fig = Figure(1,'size',[40 40]);
            
            overlap = 100*cellfun(@getOverlap, C0, S1);
            connectivity = 100*(1-cellfun(@cove.sparsity, S1));
            scatter(connectivity,overlap,17,nCells,'filled')
            colormap(flipud(hot)/1.5)
            hold on, plot(connectivity(hix2),overlap(hix2),'bo')
            set(gca,'XTick',0:20:100)
            xlabel '% connectivity'
            set(gca,'YTick',0:20:100,'XTick',0:10:100)
            ylabel '% overlap with high corrs    '
            ylim(ylim.*[0 1])
            xlim([0 max(connectivity)])
            fprintf('mean overlap: %2.1f%%\n', mean(overlap))
            
            fig.cleanup
            fig.save(fname)
            
            
            % Panel F: % negative interactions
            fname = fullfile(covest.plots.figPath, 'Fig5-F.eps');
            fig = Figure(1,'size',[40 40]);
            
            pcentNegative = 100*cellfun(@fracNegative, S1);
            scatter(connectivity,pcentNegative,17,nCells,'filled')
            colormap(flipud(hot)/1.5)
            hold on, plot(connectivity(hix2),pcentNegative(hix2),'bo')
            set(gca,'XTick',0:10:100)
            xlabel '% connectivity'
            set(gca,'YTick',0:20:100)
            ylabel '% neg. interactions'
            ylim(ylim.*[0 1])
            xlim([0 max(connectivity)])
            fprintf('mean percent negative: %2.1f%%\n', mean(pcentNegative))
            
            fig.cleanup
            fig.save(fname)
            
            function overlap = getOverlap(C,S)
                p = size(C);
                [i,j] = ndgrid(1:p,1:p);
                thresh = quantile(abs(C(i<j)),cove.sparsity(S));
                C = C.*(abs(C)>thresh);
                overlap = sum(logical(C(i<j) & S(i<j)))/sum(logical(S(i<j)));
            end
            
            function frac = fracNegative(S)
                p = size(S);
                [i,j] = ndgrid(1:p,1:p);
                frac = sum(S(i<j)>0)/sum(logical(S(i<j)));
            end
            
        end
        
        
        function fig6_dist
            % Figure 6, panels B, C, E, F
            for vertical = [false true]  % false = panels B, D,  false = C, F
                c = covest.CovMatrix & 'nfolds=1';
                c0 = c.pro('method->m0','cov_matrix->c0') & 'm0=0';
                c1 = c.pro('method->m1','cov_matrix','sparse') & 'm1=90';
                rel = c0*c1*covest.Traces*covest.ActiveCells;
                [xyz,selection,C0,C1,S] = rel.fetchn('cell_xyz','selection','c0','cov_matrix','sparse');
                
                % remove invactive cells
                xyz = cellfun(@(xyz,selection) xyz(selection,:), xyz, selection, 'uni', false);
                
                % convert to correlations and partial correlations, respectively
                C0 = cellfun(@corrcov, C0, 'uni', false);  % sample correlations
                C1 = cellfun(@(C) -corrcov(inv(C)), C1, 'uni', false);  % regularized partial correlations
                
                % compute average correlations
                avgC0 = cellfun(@(C) mean(offDiag(C)), C0);
                avgC1 = cellfun(@(C) mean(offDiag(C)), C1);
                avgConn = cellfun(@(S) mean(logical(offDiag(S))), S);  % average connectivity
                
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
                    edges = [0 25 60  120];
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
                
                % bin distances
                nBins = length(edges)-1;
                D = cellfun(@(oriDiffs) sum(bsxfun(@ge,oriDiffs,edges),2), D, 'uni', false);
                
                % exclude sites that don't have enough cell pairs in each bin
                hasEnough = 100 < cellfun(@(D) min(hist(D,1:nBins)), D);
                fprintf('Qualifying sites n=%d\n', sum(hasEnough))
                D = D(hasEnough);
                C0 = C0(hasEnough);
                C1 = C1(hasEnough);
                S = S(hasEnough);
                avgC0 = avgC0(hasEnough);
                avgC1 = avgC1(hasEnough);
                avgConn = avgConn(hasEnough);
                clear hasEnough
                
                % panel B or C: average correlations vs distance
                if vertical
                    fig = Figure(1,'size',[40,30]);
                    fname = fullfile(covest.plots.figPath, 'Fig6-C.eps');
                else
                    fig = Figure(1,'size',[45,30]);
                    fname = fullfile(covest.plots.figPath, 'Fig6-B.eps');
                end
                
                C0 = cellfun(@(C,D) accumarray(D,C,[nBins 1],@mean), C0, D, 'uni',false);
                C1 = cellfun(@(C,D) accumarray(D,C,[nBins 1],@mean), C1, D, 'uni',false);
                C0 = bsxfun(@rdivide, cat(2,C0{:})', avgC0);
                C1 = bsxfun(@rdivide, cat(2,C1{:})', avgC1);
                clear avgC0 avgC1
                
                ticks = conv(edges,[.5 .5],'valid')';
                
                errorbar(ticks,mean(C0),std(C0)/sqrt(size(C0,1)), 'k.-')
                hold on
                errorbar(ticks,mean(C1),std(C1)/sqrt(size(C1,1)), 'r.-')
                plot(ticks,mean(C0), 'k.-', 'LineWidth', 1)
                plot(ticks,mean(C1), 'r.-', 'LineWidth', 1)
                plot(xlim,[1 1],'k:')
                hold off
                set(gca,'XTick',edges(1:end-1),'YTick',[0 1])
                xlim([0 edges(end)])
                ylim(ylim.*[0 1])
                if vertical
                    xlabel 'depth displacement (\mum)  '
                else
                    xlabel 'lateral displacement (\mum)  '
                end
                ylabel 'avg. norm. corr.'
                
                fig.cleanup
                fig.save(fname)
                
                
                % select only sites with sufficient connectivity
                clear C1 C0
                isConn = avgConn > 0.05;
                fprintf('Qualifying sites for connectivity vs ori: n=%d\n', sum(isConn))
                D = D(isConn);
                S = S(isConn);
                avgConn = avgConn(isConn);
                
                % panels E and F: connectivity vs ori tuning
                if vertical
                    fig = Figure(1,'size',[40,30]);
                    fname = fullfile(covest.plots.figPath, 'Fig6-F.eps');
                else
                    fig = Figure(1,'size',[45,30]);
                    fname = fullfile(covest.plots.figPath, 'Fig6-E.eps');
                end
                posColor = [0 .5 0];
                negColor = [.5 0 0];
                
                neg = cellfun(@(S,D) accumarray(D,S>0,[nBins 1], @mean), S, D, 'uni',false);
                pos = cellfun(@(S,D) accumarray(D,S<0,[nBins 1], @mean), S, D, 'uni',false);
                neg = bsxfun(@rdivide, cat(2,neg{:})', avgConn);
                pos = bsxfun(@rdivide, cat(2,pos{:})', avgConn);
                clear avgConn
                
                errorbar(ticks,mean(pos),std(pos)/sqrt(size(pos,1)), '.-', 'Color',posColor)
                hold on
                errorbar(ticks,mean(neg),std(neg)/sqrt(size(neg,1)), '.-', 'Color',negColor)
                
                plot(ticks,mean(pos), '.-', 'LineWidth', 1, 'Color', posColor)
                plot(ticks,mean(neg), '.-', 'LineWidth', 1, 'Color', negColor)
                plot(xlim,[1 1],'k:')
                hold off
                set(gca,'XTick',edges(1:end-1),'YTick',[0 1])
                xlim([0 edges(end)])
                ylim(ylim.*[0 1])
                if vertical
                    xlabel 'depth displacement (\mum)'
                else
                    xlabel 'lateral displacement (\mum)'
                end
                ylabel 'rel. connectivity'
                
                fig.cleanup
                fig.save(fname)
            end
            
            
            function D = distMatrix(xyz)
                p = size(xyz,1);
                [i,j] = ndgrid(1:p,1:p);
                D = reshape(sqrt(sum((xyz(i,:)-xyz(j,:)).^2,2)),p,p);
            end
        end
        
        
        function fig6_ori
            % Figure 6, panels A and D.
            
            % get all data
            c = covest.CovMatrix & 'nfolds=1';
            c0 = c.pro('method->m0', 'cov_matrix->c0') & 'm0=0';
            c1 = c.pro('method->m1','sparse','cov_matrix') & 'm1=90';
            rel = c0*c1*covest.ActiveCells*pro(covest.OriTuning,'high_repeats->hrepeats','*');
            [pval,ori,selection,C0,C1,S] = rel.fetchn('von_p_value','von_pref','selection','c0','cov_matrix','sparse');
            
            % remove inactive cells
            ori = cellfun(@(ori,selection) ori(selection)*180/pi, ori, selection,'uni',false);
            pval = cellfun(@(pval,selection) pval(selection), pval, selection,'uni',false);
            clear selection
            
            % convert to correlations and partial correlations, respectively
            C0 = cellfun(@corrcov, C0, 'uni', false);  % sample correlations
            C1 = cellfun(@(C) -corrcov(inv(C)), C1, 'uni', false);  % regularized partial correlations
            
            % compute average correlations
            avgC0 = cellfun(@(C) mean(offDiag(C)), C0);
            avgC1 = cellfun(@(C) mean(offDiag(C)), C1);
            avgConn = cellfun(@(S) mean(logical(offDiag(S))), S);  % average connectivity
            
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
            
            % exclude sites that don't have enough tuned cell pairs in each bin
            hasEnough = 100 < cellfun(@(D) min(hist(D,1:nBins)), D);
            fprintf('Qualifying sites n=%d\n', sum(hasEnough))
            D = D(hasEnough);
            C0 = C0(hasEnough);
            C1 = C1(hasEnough);
            S = S(hasEnough);
            avgC0 = avgC0(hasEnough);
            avgC1 = avgC1(hasEnough);
            avgConn = avgConn(hasEnough);
            clear hasEnough
            
            % panel A: average correlations vs ori tuning
            fig = Figure(1,'size',[40,30]);
            fname = fullfile(covest.plots.figPath, 'Fig6-A.eps');
            
            C0 = cellfun(@(C,D) accumarray(D,C,[nBins 1],@mean), C0, D, 'uni',false);
            C1 = cellfun(@(C,D) accumarray(D,C,[nBins 1],@mean), C1, D, 'uni',false);
            C0 = bsxfun(@rdivide, cat(2,C0{:})', avgC0);
            C1 = bsxfun(@rdivide, cat(2,C1{:})', avgC1);
            clear avgC0 avgC1
            
            ticks = conv(edges,[.5 .5],'valid')';
            
            errorbar(ticks,mean(C0),std(C0)/sqrt(size(C0,1)), 'k.-')
            hold on
            errorbar(ticks,mean(C1),std(C1)/sqrt(size(C1,1)), 'r.-')
            plot(ticks,mean(C0), 'k.-', 'LineWidth', 1)
            plot(ticks,mean(C1), 'r.-', 'LineWidth', 1)
            plot(xlim,[1 1],'k:')
            hold off
            set(gca,'XTick',edges,'YTick',[0 1])
            xlim([0 90])
            ylim(ylim.*[0 1])
            xlabel '\Deltaori (degrees)'
            ylabel 'avg. norm. corr.'
            fig.cleanup
            fig.save(fname)
            
            % select only sites with sufficient connectivity
            isConn = avgConn > 0.05;
            fprintf('Qualifying sites for connectivity vs ori: n=%d\n', sum(isConn))
            D = D(isConn);
            S = S(isConn);
            avgConn = avgConn(isConn);
            
            % panel D: connectivity vs ori tuning
            fig = Figure(1,'size',[40,30]);
            fname = fullfile(covest.plots.figPath, 'Fig6-D.eps');
            posColor = [0 .5 0];
            negColor = [.5 0 0];
            
            neg = cellfun(@(S,D) accumarray(D,S>0,[nBins 1], @mean), S, D, 'uni',false);
            pos = cellfun(@(S,D) accumarray(D,S<0,[nBins 1], @mean), S, D, 'uni',false);
            neg = bsxfun(@rdivide, cat(2,neg{:})', avgConn);
            pos = bsxfun(@rdivide, cat(2,pos{:})', avgConn);
            clear avgConn
            
            errorbar(ticks,mean(pos),std(pos)/sqrt(size(pos,1)), '.-', 'Color',posColor)
            hold on
            errorbar(ticks,mean(neg),std(neg)/sqrt(size(neg,1)), '.-', 'Color',negColor)
            
            plot(ticks,mean(pos), '.-', 'LineWidth', 1, 'Color', posColor)
            plot(ticks,mean(neg), '.-', 'LineWidth', 1, 'Color', negColor)
            plot(xlim,[1 1],'k:')
            hold off
            set(gca,'XTick',edges,'YTick',[0 1])
            xlim([0 90])
            ylim(ylim.*[0 1])
            xlabel '\Deltaori (degrees)'
            ylabel 'rel. connectivity'
            
            fig.cleanup
            fig.save(fname)
        end
        
        function supp1A
            alphas = [.05 .2 .3 .4 .5 .6 .8 1];
            L = zeros(0,length(alphas));
            nfolds = 10;
            for key = fetch(covest.ActiveCells & 'preprocess_method_num=5' & 'high_repeats' & 'ncells>100')'
                fprintf .
                [X,evokedBins,sel] = fetch1(covest.ActiveCells*covest.Traces & key, ...
                    'trace_segments','evoked_bins', 'selection');
                X = X(1:min(end,evokedBins),:,:,sel);
                [nBins,nCond,nTrials,nCells] = size(X); %#ok<ASGLU>
                L(end+1,:) = arrayfun(@(alpha) ...
                    mean(arrayfun(@(k) getLoss(X,k,nfolds,true,alpha), 1:nfolds)) - ...
                    mean(arrayfun(@(k) getLoss(X,k,nfolds,false), 1:nfolds)),...
                    alphas);     %#ok<AGROW>
            end
            
            fig = Figure(1, 'size', [163 35]);
            h = boxplot(L,'jitter',0,'colors','k',...
                'labels',arrayfun(@(a) sprintf('a = %1.2f',a), alphas,'uni',false),...
                'orientation','horizontal','outliersize',3);
            set(h(1:2,:),'LineStyle','-','LineWidth',.25)
            set(h(7,:),'MarkerEdgeColor','k')
            xlabel nats/cell/bin
            set(gca,'YColor',[1 1 1],'YDir','reverse')
            hold on
            plot([0 0],ylim,'k:')
            hold off
            set(gca,'Position',[0.1 0.3 0.89 0.7])
            
            fig.cleanup
            fig.save(fullfile(covest.plots.figPath,'Supp1A.eps'))
            
            
            
            
            function L = getLoss(X, k, nfolds, binwise, alpha)
                
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
                
                L = loss(C, cove.cov(reshape(Y,[],nCells))) + nanmean(log(V(:)));
            end
        end
        
        
        
        function supp1BC
            % Supplementary Figure 1, panels B and C
            % assess effect of common variance assumption on limited-range correlations
            
            C1 = zeros(0,2);  % average correlation of similarly tuned cells
            C2 = zeros(0,2);  % average correlation of differently tuned cells
            oriRel = pro(covest.OriTuning,'high_repeats->temp','*');
            for key = fetch(covest.ActiveCells*oriRel & 'ncells>100')'
                [X,evokedBins,sel,ori,pval] = fetch1(covest.ActiveCells*covest.Traces*oriRel & key, ...
                    'trace_segments','evoked_bins', 'selection','von_pref','von_p_value');
                X = double(X(1:min(end,evokedBins),:,:,sel));  % only evoked bins
                ori = ori(sel)*180/pi;
                pval = pval(sel);
                c0 = corrVar(X,'bin');
                c1 = corrVar(X,'common');
                
                q = [nan 1];
                ori(pval>0.05)=nan;  % exclude untuned cells
                od = oriDiff(ori,ori');
                
                % average correlations between differently tuned cells
                ix = q((~isnan(od) & od < 15)+1);
                [cc1,n] = cellfun(@(c) avgOffDiag(c.*ix), {c0 c1});
                if any(n<300), continue, end
                
                % average correlations between similarly tuned cells
                ix = q((~isnan(od) & od > 45)+1);
                [cc2,n] = cellfun(@(c) avgOffDiag(c.*ix), {c0 c1});
                if any(n<300), continue, end
                
                C1(end+1,:) = cc1; %#ok<AGROW>
                C2(end+1,:) = cc2; %#ok<AGROW>
            end
            fig = Figure(1,'size',[80 25]);
            h = boxplot([C1(:,1) C2(:,1)],'jitter',0,'colors','k',...
                'labels',{'<15','>45'},'orientation','horizontal','outliersize',3);
            set(h(1:2,:),'LineStyle','-','LineWidth',.25)
            set(h(7,:),'MarkerEdgeColor','k')
            xlabel 'avg corr'
            xlim(xlim.*[0 1])
            fig.cleanup
            fig.save(fullfile(covest.plots.figPath,'Supp1B.eps'))
            
            fig = Figure(1,'size',[80 25]);
            h = boxplot(C1-C2,'jitter',0,'colors','k',...
                'labels',{'per bin' 'common'},'orientation','horizontal','outliersize',3);
            set(h(1:2,:),'LineStyle','-','LineWidth',.25)
            set(h(7,:),'MarkerEdgeColor','k')
            xlim(xlim.*[0 1])
            xlabel '\Delta avg corr'
            fig.cleanup
            fig.save(fullfile(covest.plots.figPath,'Supp1C.eps'))            
            
        end
    end
end


function s=nozero(f)
% remove leading zeros in decimal fractions
if isscalar(f)
    s = sprintf('%g',f);
    if strncmp(s,'-0.',3)
        s(2)='';
    elseif strncmp(s,'0.',2)
        s(1)='';
    end
else
    s = arrayfun(@(f) nozero(f), f, 'uni', false);
end
end


function ret = oriDiffMatrix(ori)
% make a matrix of orientation differences
ret = oriDiff(ori,ori');
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

function [ret,n] = avgOffDiag(C)
ret = nanmean(offDiag(C));
if nargout>1
    n = sum(~isnan(offDiag(C)));
end
end