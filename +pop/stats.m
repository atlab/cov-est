classdef stats
    
    methods(Static)
        
        function overlap
            %histogram of overlap with thresholded correlations
            interest = 4328;
            r = pop.AodCovEstimate2;
            r1 = pro(r & 'cov_estim_num=0','cov_matrix->c0','cov_estim_num->cn');
            [s,it,C0] = fetchn(r*r1 & 'cov_estim_num=14' & 'bin_opt=0', ...
                'sparse',sprintf('(mod(aod_scan_start_time,1e4)=%d)->it', interest),'c0');
            c = 1-cellfun(@covest.sparsity, s);
            nneurons = cellfun(@(C) size(C,1),C0);
            overlap = cellfun(@fracOverlap, s, C0);
            
            fig = Figure(1,'size',[40 40]);
            it = find(it);
            c = c*100;
            overlap = overlap*100;
            scatter(c,overlap,17,nneurons,'filled')
            hold on
            scatter(c(it),overlap(it),'bo')
            hold off
            colormap((flipud(hot))/1.5)
            xlabel '% connectivity'
            ylabel '% overlap with high corrs'
            ylim([0 80])
            xlim([0 28])
            set(gca, 'XTick', [0 10 20], 'YTick', [0 25 50 75])
            fig.cleanup
            fig.save('~/cov/figures/src/stats-overlap.eps')
            
            function overlap = fracOverlap(S,C0)
                thresh = covest.sparsity(S);
                p = size(C0,1);
                [i,j] = meshgrid(1:p,1:p);
                thresh = quantile(abs(C0(i>j)), thresh);
                C0 = abs(C0)>thresh;
                overlap = (1-covest.sparsity(C0.*S))/(1-covest.sparsity(S));                
            end
        end
        
        function negatives
            % histogram of the average number fractions of negative
            interest = 4328;
            [s,p,it] = fetchn(pop.AodCovEstimate2*pop.AodBinnedTraces2 & 'cov_estim_num=14' & 'bin_opt=0', ...
                'sparse','nneurons',sprintf('(mod(aod_scan_start_time,1e4)=%d)->it', interest));
            f = cellfun(@fracNegative, s);
            c = 1-cellfun(@covest.sparsity, s);
            c = c*100;
            f = f*100;
            
            fig = Figure(1,'size',[40 40]);
            it = find(it);
            scatter(c,f,17,p,'filled')
            hold on
            scatter(c(it),f(it),'bo')
            hold off
            colormap((flipud(hot))/1.5)
            xlabel '% connectivity'
            ylabel '% negative'
            ylim([0 50])
            xlim([0 28])
            set(gca, 'XTick', [0 10 20], 'YTick', [0 25 50])
            fig.cleanup

            fig.save('~/cov/figures/src/stats-frac-negative.eps')
            
            function f = fracNegative(s)
                f = nan;
                if sum(~~s(:))>1000
                    f = (1-covest.sparsity(s>0))/(1-covest.sparsity(s));
                end
            end
            
        end
        
        function latent
            % number of latent variables vs distance from center
            keys = fetch(pop.AodCovEstimate2 & 'cov_estim_num=14' & 'bin_opt=0');
            s = arrayfun(@getInfo, keys);
            scanNum = 4328;
            interest = find(mod([keys.aod_scan_start_time],1e4)==scanNum);
            fig = Figure(1,'size',[40 40]);
            scatter(abs(s(interest).distanceFromCenter),s(interest).latentWeight,10,'k.')
            refline
            xlim([0 150])
            fig.cleanup
            fig.save(sprintf('~/cov/figures/src/stats-latentweight-vs-position-%d.eps',scanNum))
           
            function ret = getInfo(key)
                cellnums= fetch1(pop.AodBinnedTraces2 & key, 'cellnums');
                [x,y,z,xyzcellnums] = fetchn(aod.Traces*aod.UniqueCell & key ...
                    ,'x','y','z','cell_num');
                [S,L] = fetch1(pop.AodCovEstimate2 & key, 'sparse','lowrank');
                P0 = S-L*L';
                d = sqrt(diag(diag(P0)));
                ix = ismember(xyzcellnums,cellnums);
                xyzcellnums = xyzcellnums(ix);
                assert(isequal(xyzcellnums,cellnums(:)))
                x = x(ix);
                y = y(ix);
                z = z(ix);
                center = mean([x y z]);
                p = size(S,1);
                ret.nodeDegree = sum((S~=0).*(1-eye(p)));
                ret.depth = z - center(3);
                ret.lateralDistance = sqrt(sum(bsxfun(@minus,[x y],center(1:2)).^2,2));
                ret.distanceFromCenter = sqrt(sum(bsxfun(@minus, [x y z],center).^2,2));  
                ret.latentWeight = diag(d\L*L'/d);
                
            end
        end
        
        
        function sparseLatent
            [Sigma0,k1] = fetchn(pop.AodCovEstimate2 & 'cov_estim_num=0','cov_matrix');
            [C,S,L,hypers,nneurons,sparsity,stime,k2] = fetchn(pop.AodCovEstimate2*pop.AodBinnedTraces2 & 'cov_estim_num=14', ...
                'cov_matrix','sparse','lowrank','hypers','nneurons','sparsity','aod_scan_start_time');
            assert(isequal(rmfield(k1,'cov_estim_num'),rmfield(k2,'cov_estim_num')))
            nlatent = cellfun(@(hypers) hypers(1), hypers);
            avgCorr = cellfun(@(s) avgOffDiag(corrcov(s)), Sigma0);
            avgPCorr = cellfun(@(s,l) avgOffDiag(-corrcov(s-l*l')), S, L); 
            degree = (nneurons-1).*(1-sparsity);
            
            interest = 4328;
            hix = find(mod(stime,1e4) == interest);            
            close all
            
            fig = Figure(1,'size',[40 40]);
            scatter(avgCorr,avgPCorr,17,nneurons,'o','filled');
            colormap((flipud(hot))/1.5)
            hold on
            plot(avgCorr(hix),avgPCorr(hix),'bo')
            xlim([0 .05])
            ylim([0 .005])
            xticks = [0 .03];
            yticks = xticks/10;
            set(gca,'XTick',xticks,'XTickLabel',nozero(xticks),...
                'YTick',yticks,'YTickLabel',nozero(yticks))
            xlabel 'avg sample corr'
            ylabel 'avg reg partial corr' 
            fig.cleanup
            fig.save('~/cov/figures/src/stats-corrs-vs-pcorrs.eps')
            
            
            fig = Figure(1,'size',[40 40]);
            scatter(nlatent,degree,17,nneurons,'filled')
            colormap((flipud(hot))/1.5)
            hold on, plot(nlatent(hix),degree(hix),'bo')
            set(gca,'YTick',[0 30 60],'XTick',[0 20 40])
            xlabel '# latent'
            ylabel 'avg node degree'
            fig.cleanup
            fig.save('~/cov/figures/src/stats-latent-vs-ndegree.eps')
            
            fig = Figure(1,'size',[40 40]);
            scatter(nneurons,nlatent,17,nneurons,'filled')
            colormap((flipud(hot))/1.5)
            hold on, plot(nneurons(hix),nlatent(hix),'bo')
            set(gca,'XTick',[0 150 300])
            xlabel '# neurons'
            set(gca,'YTick',[0 20 40])
            ylabel '# latent'
            fig.cleanup
            fig.save('~/cov/figures/src/stats-nneurons-vs-nlatent.eps')
            
            fig = Figure(1,'size',[40 40]);
            scatter(nneurons,degree,17,nneurons,'filled')
            colormap((flipud(hot))/1.5)
            hold on, plot(nneurons(hix),degree(hix),'bo')
            set(gca,'XTick',[0 150 300])
            xlabel '# neurons'
            set(gca,'YTick',[0 30 60])
            ylabel 'avg node degree'
            fig.cleanup
            fig.save('~/cov/figures/src/stats-nneurons-vs-ndegree.eps')
            
            
            function avg = avgOffDiag(C)
                p = size(C,1);
                [i,j] = ndgrid(1:p,1:p);
                avg = mean(C(i>j));
            end

        end
        
        
        function ori
            scanNums = fetchn(pop.AodCovEstimate2 ...
                & 'cov_estim_num=14','mod(aod_scan_start_time,10000)->sn');
            oriBounds = [0 15 45 90];
            interest = 4328;
            latBounds = [0 25 75 150 250];
            depBounds = [0 25 60 100];
            
            accum = [];
            for sn = scanNums'
                [avgCorr,avgTunedCorr,avgPCorr,avgTunedPCorr,numTuned,numPos,numNeg,key] = pop.ori(sn, oriBounds);
                fracConnected = 1-fetch1(pop.AodCovEstimate2 & key, 'sparsity');
                fracPosConnected = 1-covest.sparsity(fetch1(pop.AodCovEstimate2 & key,'sparse')<0);
                [depCorr,depPCorr,depTotal,depPos,depNeg] = pop.dist(sn,depBounds,[0 25]);
                [latCorr,latPCorr,latTotal,latPos,latNeg] = pop.dist(sn,[0 25],latBounds);
                accum = [accum struct(...
                    'interest', sn==interest, ...
                    'avgCorr', avgCorr, ...
                    'avgTunedCorr', avgTunedCorr', ...
                    'avgPCorr',avgPCorr, ...
                    'avgTunedPCorr', avgTunedPCorr', ...
                    'numTuned', numTuned', ...
                    'numPos', numPos', ...
                    'numNeg', numNeg', ...
                    'latCorr', latCorr(:), ...
                    'latPCorr', latPCorr(:), ...
                    'depCorr', depCorr(:), ...
                    'depPCorr', depPCorr(:), ...
                    'depTotal', depTotal(:), ...
                    'depPos', depPos(:), ...
                    'depNeg', depNeg(:), ...
                    'latTotal', latTotal(:), ...
                    'latPos', latPos(:), ...
                    'latNeg', latNeg(:), ...
                    'fracConnected', fracConnected, ...
                    'fracPosConnected', fracPosConnected...
                    )]; %#ok<AGROW>         
                fprintf .
            end
            fprintf \n
            
            
            %%%%%%%%%%%%%%%%%%% ORIENTATION %%%%%%%%%%%%%%%%%%%%%
            
            % ori - corr
            fig = Figure(1,'size',[35,30]);
            ix = [accum.interest];
            ticks = conv(oriBounds,[.5 .5],'valid')';            
            plot(ticks,accum(ix).avgTunedCorr/accum(ix).avgCorr,'k.-','LineWidth',1)
            set(gca,'XTick',oriBounds,'YTick',[0 1])
            hold on
            plot(ticks,accum(ix).avgTunedPCorr/accum(ix).avgPCorr,'r.-','LineWidth',1);
            plot([0 90],[1 1],'k:')
            hold off
            ylabel 'avg norm corr'
            xlabel '\Deltaori (degrees)'
            xlim([0 90])
            ylim([0 4.2])
            fig.cleanup
            fig.save(sprintf('~/cov/figures/src/summary-corr-vs-ori-%4d.eps',interest))            
            
            fig = Figure(1,'size',[35 30]);
            ix = min([accum.numTuned])>10;
            c = [accum(ix).avgTunedCorr]/diag([accum(ix).avgCorr]);
            errorbar(ticks,mean(c,2),std(c,[],2)/sqrt(sum(ix)),'k.-');
            hold on
            plot(ticks,mean(c,2),'k.-','LineWidth',1);
            c = [accum(ix).avgTunedPCorr]/diag([accum(ix).avgPCorr]);
            errorbar(ticks,mean(c,2),std(c,[],2)/sqrt(sum(ix)),'r.-');
            plot(ticks,mean(c,2),'r.-','LineWidth',1);
            plot(xlim,[1 1],'k:')
            hold off
            set(gca,'XTick',oriBounds)
            ylabel 'avg norm corr'
            xlabel '\Deltaori (degrees)'
            xlim([0 90])
            ylim([0 3.5])
            set(gca,'YTick',[0 1])
            fig.cleanup
            fig.save('~/cov/figures/src/summary-corr-vs-ori.eps')
            fprintf('corr vs ori  n=%d\n', sum(ix))
            
            % ori - conn
            fig = Figure(1,'size',[35,30]);
            posColor = [0 .5 0];
            negColor = [.5 0 0];
            ix = [accum.interest];
            plot(ticks,accum(ix).numPos./accum(ix).numTuned/accum(ix).fracConnected,'.-','LineWidth',1,'Color',posColor)
            hold on
            plot(ticks,accum(ix).numNeg./accum(ix).numTuned/accum(ix).fracConnected,'.-','LineWidth',1,'Color',negColor)
            plot([0 90],[1 1],'k:')
            set(gca,'XTick',oriBounds,'YTick',[0 1])
            hold off
            ylabel 'rel connectivity'
            xlabel '\Deltaori (degrees)'
            xlim([0 90])
            ylim([0 2.0])
            fig.cleanup
            fig.save(sprintf('~/cov/figures/src/summary-conn-vs-ori-%4d.eps',interest))     
            
            
            fig = Figure(1,'size',[35 30]);
            ix = max([accum.numNeg]+[accum.numPos])>20 & min([accum.numTuned])>50;
            c = ([accum(ix).numPos]./[accum(ix).numTuned])/diag([accum(ix).fracConnected]);
            errorbar(ticks,mean(c,2),std(c,[],2)/sqrt(sum(ix)),'.-','Color',posColor)
            hold on
            plot(ticks,mean(c,2),'.-','Color',posColor,'LineWidth',1);
            c = ([accum(ix).numNeg]./[accum(ix).numTuned])/diag([accum(ix).fracConnected]);
            errorbar(ticks,mean(c,2),std(c,[],2)/sqrt(sum(ix)),'Color',negColor)
            plot(ticks,mean(c,2),'Color',negColor,'LineWidth',1)
            plot(xlim,[1 1],'k:')
            hold off
            set(gca,'XTick',oriBounds)
            ylabel 'rel connectivity'
            xlabel '\Deltaori (degrees)'
            xlim([0 90])
            ylim([0 2.0])
            set(gca,'YTick',[0 1])
            fig.cleanup
            fig.save('~/cov/figures/src/summary-conn-vs-ori.eps')
            fprintf('conn vs ori  n=%d\n', sum(ix))
            

            %%%%%%%%%%%%%%%%%%%% DEPTH %%%%%%%%%%%%%%%%%%%%%%%%
            
            % depth - corr
            fig = Figure(1,'size',[35 30]);
            ticks = conv(depBounds,[.5 .5],'valid');            
            ix = [accum.interest];
            plot(ticks,accum(ix).depCorr/accum(ix).avgCorr,'k.-','LineWidth',1)
            hold on
            plot(ticks,accum(ix).depPCorr/accum(ix).avgPCorr,'r.-','LineWidth',1)
            plot(xlim,[1 1],'k:')
            hold off
            set(gca,'XTick',depBounds(1:end-1),'YTick',[0 1])
            xlim([0 ticks(end)+20])
            ylim([0 4.2])      
            xlabel 'depth displacement (\mum)    '
            ylabel 'avg norm corr'
            fig.cleanup
            fig.save(sprintf('~/cov/figures/src/summary-corr-vs-depth-%04d.eps',interest))
            
            
            fig = Figure(1,'size',[35 30]);
            ix = min([accum.depTotal])>12;
            c = [accum(ix).depCorr]/diag([accum(ix).avgCorr]);
            errorbar(ticks,mean(c,2),std(c,[],2)/sqrt(sum(ix)),'k.-')
            hold on
            plot(ticks,mean(c,2),'k.-','LineWidth',1)
            errorbar(ticks,mean(c,2),std(c,[],2)/sqrt(sum(ix)),'k.-')
            c = [accum(ix).depPCorr]/diag([accum(ix).avgPCorr]);
            errorbar(ticks,mean(c,2),std(c,[],2)/sqrt(sum(ix)),'r.-')
            plot(ticks,mean(c,2),'r.-','LineWidth',1)
            plot(xlim,[1 1],'k:')
            hold off
            set(gca,'XTick',depBounds(1:end-1),'YTick',[0 1])
            xlim([0 ticks(end)+20])
            ylim([0 3.5])      
            xlabel 'depth displacement (\mum)    '
            ylabel 'avg norm corr'
            fig.cleanup
            fig.save('~/cov/figures/src/summary-corr-vs-depth.eps')
            fprintf('corr vs depth  n=%d\n', sum(ix))
            
            
            % depth - conn
            fig = Figure(1,'size',[35 30]);
            ticks = conv(depBounds,[.5 .5],'valid');            
            ix = [accum.interest];
            plot(ticks,(accum(ix).depPos./accum(ix).depTotal)/accum(ix).fracConnected,...
                '.-','LineWidth',1,'Color',posColor)
            hold on
            plot(ticks,(accum(ix).depNeg./accum(ix).depTotal)/accum(ix).fracConnected,...
                '.-','LineWidth',1,'Color',negColor)
            plot(xlim,[1 1],'k:')
            hold off
            set(gca,'XTick',depBounds(1:end-1),'YTick',[0 1])
            xlim([0 ticks(end)+20])
            ylim([0 2.0])      
            xlabel 'depth displacement (\mum)    '
            ylabel 'rel connectivity'
            fig.cleanup
            fig.save(sprintf('~/cov/figures/src/summary-conn-vs-depth-%04d.eps',interest))
            
            
            fig = Figure(1,'size',[35 30]);
            ix = max([accum.depNeg]+[accum.depPos])>20 & min([accum.depTotal])>50;
            c = ([accum(ix).depPos]./[accum(ix).depTotal])/diag([accum(ix).fracConnected]);
            errorbar(ticks,mean(c,2),std(c,[],2)/sqrt(sum(ix)),'.-','Color',posColor)
            hold on
            plot(ticks,mean(c,2),'.-','LineWidth',1,'color',posColor)
            c = ([accum(ix).depNeg]./[accum(ix).depTotal])/diag([accum(ix).fracConnected]);
            errorbar(ticks,mean(c,2),std(c,[],2)/sqrt(sum(ix)),'.-','Color',negColor)
            plot(ticks,mean(c,2),'.-','LineWidth',1,'Color',negColor)
            plot(xlim,[1 1],'k:')
            hold off
            set(gca,'XTick',depBounds(1:end-1),'YTick',[0 1])
            xlim([0 ticks(end)+20])
            ylim([0 2.0])      
            xlabel 'depth displacement (\mum)    '
            ylabel 'rel connectivity'
            fig.cleanup
            fig.save('~/cov/figures/src/summary-conn-vs-depth.eps')
            fprintf('conn vs depth  n=%d\n', sum(ix))

            %%%%%%%%%%%%%%%%%% LATERAL %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % lateral - ori
            fig = Figure(1,'size',[40 30]);
            ticks = conv(latBounds,[.5 .5],'valid');
            ix = [accum.interest];
            plot(ticks,accum(ix).latCorr/accum(ix).avgCorr,'k.-','LineWidth',1)
            hold on
            plot(ticks,accum(ix).latPCorr/accum(ix).avgPCorr,'r.-','LineWidth',1)
            plot(xlim,[1 1],'k:')
            hold off
            set(gca,'XTick',latBounds(1:end-1),'YTick',[0 1])
            xlim([0 ticks(end)+20])
            ylim([0 4.2])      
            xlabel 'lateral displacement (\mum)    '
            ylabel 'avg norm corr'
            fig.cleanup
            fig.save(sprintf('~/cov/figures/src/summary-corr-vs-lat-%04d.eps',interest))
            
            
            fig = Figure(1,'size',[40 30]);
            ix = min([accum.latTotal])>12;
            c = [accum(ix).latCorr]/diag([accum(ix).avgCorr]);
            errorbar(ticks,mean(c,2),std(c,[],2)/sqrt(sum(ix)),'k.-')
            hold on
            plot(ticks,mean(c,2),'k.-','LineWidth',1)
            c = [accum(ix).latPCorr]/diag([accum(ix).avgPCorr]);
            errorbar(ticks,mean(c,2),std(c,[],2)/sqrt(sum(ix)),'r.-')
            plot(ticks,mean(c,2),'r.-','LineWidth',1)
            plot(xlim,[1 1],'k:')
            hold off
            set(gca,'XTick',latBounds(1:end-1),'YTick',[0 1])
            xlim([0 ticks(end)+20])
            ylim([0 3.5])      
            xlabel 'lateral displacement (\mum)    '
            ylabel 'avg norm corr'
            fig.cleanup
            fig.save('~/cov/figures/src/summary-corr-vs-lat.eps')
            fprintf('corr vs lat  n=%d\n', sum(ix))
            
                        %%%% LEGEND %%%%
            fig = Figure(1,'size',[30 30]);
            plot([-1 0 1],1*[1 1 1],'k.-','LineWidth',1)
            hold on
            plot([-1 0 1], -1*[1 1 1],'r.-','LineWidth',1)
            hold off
            axis off
            set(gca,'Position',[.3 .4 .4 .2])

            fig.cleanup
            fig.save('~/cov/figures/src/summary-corr-legend.eps')


            
            
            
            % lateral - conn
            fig = Figure(1,'size',[40 30]);
            ticks = conv(latBounds,[.5 .5],'valid');
            ix = [accum.interest];
            plot(ticks,(accum(ix).latPos./accum(ix).latTotal)/accum(ix).fracConnected,...
                '.-','LineWidth',1,'Color',posColor)
            hold on
            plot(ticks,(accum(ix).latNeg./accum(ix).latTotal)/accum(ix).fracConnected,...
                '.-','LineWidth',1,'Color',negColor)
            plot(xlim,[1 1],'k:')
            hold off
            set(gca,'XTick',latBounds(1:end-1),'YTick',[0 1])
            xlim([0 ticks(end)+20])
            ylim([0 2.0])
            xlabel 'lateral displacement (\mum)    '
            ylabel 'rel connectivity'
            fig.cleanup
            fig.save(sprintf('~/cov/figures/src/summary-conn-vs-lat-%04d.eps',interest))
            
            
            fig = Figure(1,'size',[40 30]);
            ix = max([accum.latNeg]+[accum.latPos])>20 & min([accum.latTotal])>50;
            c = ([accum(ix).latPos]./[accum(ix).latTotal])/diag([accum(ix).fracConnected]);
            errorbar(ticks,mean(c,2),std(c,[],2)/sqrt(sum(ix)),'.-','Color',posColor)
            hold on
            plot(ticks,mean(c,2),'.-','LineWidth',1,'color',posColor)
            c = ([accum(ix).latNeg]./[accum(ix).latTotal])/diag([accum(ix).fracConnected]);
            errorbar(ticks,mean(c,2),std(c,[],2)/sqrt(sum(ix)),'.-','Color',negColor)
            plot(ticks,mean(c,2),'.-','LineWidth',1,'Color',negColor)
            plot(xlim,[1 1],'k:')
            hold off
            set(gca,'XTick',latBounds(1:end-1),'YTick',[0 1])
            xlim([0 ticks(end)+20])
            ylim([0 2.0])      
            xlabel 'lateral displacement (\mum)    '
            ylabel 'rel connectivity'
            fig.cleanup
            fig.save('~/cov/figures/src/summary-conn-vs-lat.eps')
            fprintf('conn vs lat  n=%d\n', sum(ix))
            
            
            %%%% LEGEND %%%%
            fig = Figure(1,'size',[30 30]);
            plot([-1 0 1],1*[1 1 1],'.-','LineWidth',1,'Color',posColor)
            hold on
            plot([-1 0 1], -1*[1 1 1],'.-','LineWidth',1,'Color',negColor)
            hold off
            axis off
            set(gca,'Position',[.3 .4 .4 .2])

            fig.cleanup
            fig.save('~/cov/figures/src/summary-conn-legend.eps')
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
