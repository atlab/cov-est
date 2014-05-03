% plots and figures for the covariance estimation paper


classdef plots
    
    properties(Constant)
        figPath = '~/cov/figures/src/'
        exampleSite = 'mod(aod_scan_start_time,10000)=4328 && high_repeats'
    end
    
    methods(Static)
        
        function fig3
            
            pairs = {
                0   90    'sample'
                10  90    'diag'
                30  90    'factor'
                80  90    'sparse'
                };
            pairs = flipud(pairs);
            
            c = covest.CovMatrix & 'nfolds>1';
            c1 = pro(c, 'method->m1','cv_loss->l1');
            c2 = pro(c, 'method->m2','cv_loss->l2');
            
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
            xlabel 'nats/cell/bin'
            ticks = 0:0.01:.1;
            set(gca,'XTick',ticks,'XTickLabel',nozero(ticks))
            set(gca,'YColor',[1 1 1])
            hold on
            plot([0 0],ylim,'k:')
            hold off
            axis tight
            set(gca,'Position',[.08 .3 0.90 0.7])
            fig.cleanup
            fig.save(fullfile(covest.plots.figPath, 'Fig3.eps'))
            
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
        end
        
        
        function fig4
            select = [covest.plots.exampleSite ' && nfolds=1'];
            C0 = fetch1(covest.CovMatrix & select & 'method=0','cov_matrix');
            [C1,S,L] = fetch1(covest.CovMatrix & select & 'method=90','cov_matrix','sparse','lowrank');
            p = size(C0,1);
            CC0 = corrcov(C0);
            
            fprintf('sparsity: %2.1f%%, avg node degree = %3.2f, low-rank=%d\n', ...
                covest.lib.sparsity(S)*100, covest.lib.nodeDegree(S), size(L,2))
            
            % correlation matrix: sample / regularized
            comboPlot(CC0,corrcov(C1),fullfile(covest.plots.figPath, 'Fig4-A.eps'))
            densityPlot(corrcov(C0),corrcov(C1),fullfile(covest.plots.figPath,'Fig4-B'))
            
            % partial correlation matrix: sample / regularized
            iC0 = -corrcov(inv(C0));
            iC0 = (~eye(p)).*iC0;
            iC = inv(C1);
            %ds = diag(sqrt(diag(iC)));
            ds = diag(sqrt(diag(S)));
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
                    threshold = quantile(abs(C1(i>j)),covest.lib.sparsity(C2));
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
                colormap(covest.lib.doppler)
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
            clf

            doFragment = false;  % when, true plot only the inner fragment
            doCorr = false; % when true, plot thresholded correlations in the fragment
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
            alphaMultiplier = 3.0*(1+doFragment);

            scanNum = fetch1(covest.Traces & covest.plots.exampleSite, 'mod(aod_scan_start_time,10000)->scannum');
            fname = fullfile(covest.plots.figPath,sprintf('network%04d',scanNum));

            if doFragment
                fragmentRadius = 48;
                paperSize = [6.5 5.0];
                lineWidth = 1;
                zm = 0.7;
                panx = -2.5;
            else
                fragmentRadius = inf;
                paperSize = [12 11];
                lineWidth = .25;
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
            C0= corrcov(fetch1(covest.CovMatrix & setfield(key, 'method', 0), 'cov_matrix')); %#ok<SFLD>
            sparsity = fetch1(covest.CovMatrix & key, 'sparsity');
            p = size(C0,1);
            [i,j] = ndgrid(1:p,1:p);
            C0 = C0.*(abs(C0)>quantile(abs(C0(i<j)), sparsity));
            
            % report everything
            fprintf('Ssparsity = %2.1f%%\n', 100*sparsity)
            fprintf('Overlap = %2.1f%%\n',  100*sum(C0(i<j) & S(i<j))/sum(~~S(i<j)))
            fprintf('Latent = %d\n', size(L,2));
            fprintf('Negative interactions = %2.1f%%\n', 100*sum(S(i<j)<0)/sum(~~S(i<j)))
            
            if doCorr
                S = C0;
                fname = fullfile(covest.plots.figPath,sprintf('network-corr%04d',scanNum));
            elseif doFragment
                fname = fullfile(covest.plots.figPath,sprintf('network-frag%04d',scanNum));
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

            
            [i,j] = ndgrid(1:size(S,1),1:size(S,2));
            
            % positive interactions
            ix = find(j(:)>i(:) & S(:)>0);
            for ix = ix'
                lineAlpha = min(1,(alphaMultiplier*1*abs(S(ix))));
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
    end
end
