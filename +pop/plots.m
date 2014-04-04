classdef plots < handle
    methods(Static)
        
        
        function ampHist(key)
            if ~nargin
                key = fetch(aod.TracePreprocessSet & 'mod(aod_scan_start_time,10000)=4328');
            end
            X = fetch1(pop.AodBinnedTraces2 & key, 'binned_traces');
            hist(X(:),-5:0.05:15)
        end
        
        
        function key = eftychios
            key.subject_id = 16;
            key.setup = 4;
            key.session_start_time  = uint64(3456412826172);
            key.aod_scan_start_time = uint64(3456417174328);
            key.preprocess_method_num = 1;
            [key.traces,keys] = fetchn(aod.TracePreprocess & aod.UniqueCell & key, 'trace');
            key.traces = [key.traces{:}];
            times = getTimes(aod.TracePreprocess & key);
            key.times = times';
            key.stim1 = getStim(key,times,true);
            key.stim2 = getStim(key,times,false);
            
            function stim = getStim(key,times, tuning)
                % limit the trial group to 2-stimulus trials
                trialGroup = stimulation.StimTrialGroup*acq.AodStimulationLink;
                if tuning
                    trialGroup = pro(trialGroup, stimulation.StimConditions, 'count(*)->ncond') & key & 'ncond>8';
                else
                    trialGroup = pro(trialGroup, stimulation.StimConditions, 'count(*)->ncond') & key & 'ncond<8';
                end
                assert(trialGroup.count == 1);
                trials = fetch(stimulation.StimTrials*trialGroup);
                conditions = fetch(stimulation.StimConditions*trialGroup,'*');
                
                % Extract segment of trials for each stimulus
                presentations = [];
                for i = 1:length(trials)
                    trial_info = fetch1(stimulation.StimTrials & trials(i), 'trial_params');
                    event = fetch(stimulation.StimTrialEvents(trials(i), 'event_type="showSubStimulus"'),'*');
                    [onset,~] = sort([event.event_time]);
                    if length(onset)==1
                        event = fetch(stimulation.StimTrialEvents(trials(i), 'event_type="clearScreen"'),'*');
                        assert(length(event)==1)
                        offset = event.event_time;
                    else
                        offset = onset + median(diff(onset));
                    end
                    
                    cond = trial_info.conditions;
                    for j=1:length(cond)
                        if offset(j) > times(1) && onset(j) < times(end);
                            s = conditions(cond(j)).condition_info;
                            s.onsetTime = onset(j);
                            s.offsetTime = offset(j);
                            presentations = [presentations; s]; %#ok<AGROW>
                        end
                    end
                end
                stim = presentations;
            end
        end
        
        
        
        
        
        function covEstimates
            
            useQuad = false;
            
            if useQuad
                fname = 'comparison-quad';
                pairs = {
                    [14 15]  .005  .004
                    [ 0 21]  .8   .5
                    [ 0 12]  .8   .5
                    [ 0 13]  .8   .5
                    [ 0 14]  .8   .5
                    [21 12]  .55  .4
                    [21 13]  .03  .02
                    [21 14]  .03  .02
                    [12 13]  .7   .5
                    [12 14]  .7   .5
                    [13 14]  .0065 .005
                    };
            else
                fname = 'comparison-norm';
                pairs = {
                    [13 15]  .005  .004
                    [14 15]  .005  .004
                    [ 0 21]  .7   .5
                    [ 0 12]  .14  .1
                    [ 0 13]  .7   .5
                    [ 0 14]  .7   .5
                    [21 12]  .7   .5
                    [21 13]  .014 .01
                    [21 14]  .014 .01
                    [12 13]  .7   .5
                    [12 14]  .7   .5
                    [13 14]  .014 .01
                    };
            end
            
            
            s = pop.CovEstimator*pop.AodBinnedTraces & 'bin_opt=0';
            if useQuad
                s1 = s.pro(pop.CVLoss, 'cov_estim_num->c1', 'avg(quadratic_loss)->l1');
                s2 = s.pro(pop.CVLoss, 'cov_estim_num->c2', 'avg(quadratic_loss)->l2');
            else
                s1 = s.pro(pop.CVLoss, 'cov_estim_num->c1', 'avg(normal_loss)->l1');
                s2 = s.pro(pop.CVLoss, 'cov_estim_num->c2', 'avg(normal_loss)->l2');
            end
            arrowKey = fetch(pop.AodBinnedTraces & ...
                'mod(aod_scan_start_time,10000)=4328' & 'bin_opt=0');
            
            nbins = 20;
            axisAlphas = [1 1];
            for i=1:size(pairs,1)
                close all
                fig = Figure(1, 'size', [40 30]);
                rr = pairs{i,2};
                bins = linspace(-rr,rr,nbins);
                r = s1*s2 & struct('c1',pairs{i,1}(1),'c2',pairs{i,1}(2));
                [l1,l2] = r.fetchn('l1','l2');
                l1(isnan(l1))=inf;
                x = l1-l2;
%                xarrow = fetch1(r & arrowKey, 'l1-l2->d');
                a = hist(min(1,x),bins);
                bar(bins,a)
                p = signrank(x);
%                 fprintf('%2d-%2d Significance %g, difference = %g, example = %g\n', ...
%                     pairs{i,1}(1),pairs{i,1}(2),p, median(x), xarrow)
                xlim([-1 1]*rr*(1+1.1/nbins))
                ylim([-.5 16])
                hold on
                plot([0 0], ylim,'k')
                plot([0 0]+median(x), ylim, ':','Color',[.6 0 0])
                hold off
                box off
                set(gca, 'YColor',[1 1 1]*0.99,'YTick',[])
                ticks = (-1:1)*pairs{i,3};
                colormap gray
                
                set(gca,'XTick',ticks,'XTickLabel',nozero(ticks))
                set(gca,'YTick',10)
                set(gca,'Position', [.05    .2    .8    .8])
                hold on
                PlotAxisAtOrigin(axisAlphas)
                %arrow([xarrow -1],[xarrow 0],'Length',10,'BaseAngle',60,'TipAngle',16);
                hold off
                
                fig.cleanup
                fig.save(sprintf('~/cov/figures/src/%s-%02d-vs-%02d.eps',fname,pairs{i,1}(1),pairs{i,1}(2)))
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
        end
        
        
        function quad
            pairs = {
                [ 0 21]  .7   .5
                [ 0 12]  .14  .1
                [ 0 13]  .7   .5
                [ 0 14]  .7   .5
                [21 12]  .7   .5
                [21 13]  .014 .01
                [21 14]  .014 .01
                [12 13]  .7   .5
                [12 14]  .7   .5
                [13 14]  .014 .01
                };
            
            s = pop.CovEstimator*pop.AodBinnedTraces2 & 'bin_opt=0';
            s1 = s.pro(pop.AodCovEstimateFold2, 'cov_estim_num->c1', 'avg(cross_loss)->l1');
            s2 = s.pro(pop.AodCovEstimateFold2, 'cov_estim_num->c2', 'avg(cross_loss)->l2');
            
            arrowKey = fetch(pop.AodBinnedTraces2 & 'mod(aod_scan_start_time,10000)=4328' & 'bin_opt=0');
            
            nbins = 20;
            axisAlphas = [1 1];
            for i=1:size(pairs,1)
                close all
                fig = Figure(1, 'size', [40 30]);
                rr = pairs{i,2};
                bins = linspace(-rr,rr,nbins);
                r = s1*s2 & struct('c1',pairs{i,1}(1),'c2',pairs{i,1}(2));
                [l1,l2] = r.fetchn('l1','l2');
                l1(isnan(l1))=inf;
                x = l1-l2;
                %xarrow = fetch1(r & arrowKey, 'l1-l2->d');
                a = hist(min(1,x),bins);
                bar(bins,a)
                p = signrank(x);
%                 fprintf('%2d-%2d Significance %g, difference = %g, example = %g\n', ...
%                     pairs{i,1}(1),pairs{i,1}(2),p, median(x), xarrow)
                xlim([-1 1]*rr*(1+1.1/nbins))
                ylim([-.5 16])
                hold on
                plot([0 0], ylim,'k')
                plot([0 0]+median(x), ylim, ':','Color',[.6 0 0])
                hold off
                box off
                set(gca, 'YColor',[1 1 1]*0.99,'YTick',[])
                ticks = (-1:1)*pairs{i,3};
                colormap gray
                
                set(gca,'XTick',ticks,'XTickLabel',nozero(ticks))
                set(gca,'YTick',10)
                set(gca,'Position', [.05    .2    .8    .8])
                hold on
                PlotAxisAtOrigin(axisAlphas)
                arrow([xarrow -1],[xarrow 0],'Length',10,'BaseAngle',60,'TipAngle',16);
                hold off
                
                fig.cleanup
                fig.save(sprintf('~/cov/figures/src/comparison-%02d-vs-%02d.eps',pairs{i,1}(1),pairs{i,1}(2)))
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
        end
        
        
        
        function aodRaster
            key = fetch(pop.AodBinnedTraces & 'mod(aod_scan_start_time,10000)=4328' & 'bin_opt=0');
            
            fig = Figure(1,'size',[45 40]);
            [X,dt] = fetch1(pop.AodBinnedTraces & key, 'binned_traces','binsize_ms');
            [n,p] = size(X);
            times = dt/1000*(0:n-1);
            mx = quantile(X(:),0.999);
            
            %plot raster - first minute
            ix = times>=0.1 & times<=60.1;
            M = X(ix,:)'/mx;
            M = max(0,(1-max(0, cat(3,M,M,M)))).^4;
            ixTraces = 251:258;
            M(ixTraces,:,[2 3])=M(ixTraces,:,[2 3])*0.8;
            image(times(ix), 1:p, M)
            set(gca,'XTick',[10 30 50], 'YTick',[1 size(X,2)])
            set(gca, 'Position', [0.16 0.12 0.82 0.84]);
            fig.cleanup
            fig.save('~/cov/figures/src/raster-first.eps')
            
            clf
            % last minute
            ix = times>=times(end)-60-1 & times<=times(end)-1;
            M = X(ix,:)'/mx;
            M = max(0,(1-max(0, cat(3,M,M,M)))).^4;
            image(times(ix), 1:p, M)
            set(gca, 'Position', [0.16 0.12 0.82 0.84]);
            set(gca,'YTick',[],'XTick',1740+[10 30 50])
            fig.cleanup
            fig.save('~/cov/figures/src/raster-last.eps')
        end
        
        
        function traces
            key = fetch(pop.AodBinnedTraces & 'mod(aod_scan_start_time,10000)=4328' & 'bin_opt=0');
            fig = Figure(1,'size',[78 40]);
            
            aodTraces = fetch(aod.TracePreprocess & setfield(key,'preprocess_method_num',3),'*');
            aodSpikes = fetch(aod.TracePreprocess & key, '*');
            cellnums = sort(fetch1(pop.AodBinnedTraces & key, 'cellnums'));
            
            aodTraces = arrayfun(@(ix) aodTraces([aodTraces.cell_num]==ix), cellnums);
            aodSpikes = arrayfun(@(ix) aodSpikes([aodSpikes.cell_num]==ix), cellnums);
            X = [aodSpikes.trace];
            mx = quantile(X(:),0.999);
            
            r = [0 60];
            traceIdx = (1:8)+250;
            step = 0.4;
            times = (0:size(aodSpikes(1).trace,1)-1)/aodSpikes(1).fs;
            ix = times>=r(1) & times<=r(2);
            X = X(ix,traceIdx)'/mx;
            X = max(0,(1-max(0, X))).^4;
            M = ones(size(X,1)*4,size(X,2),3);
            M(2:4:end,:,1) = 1;
            M(2:4:end,:,2) = X;
            M(2:4:end,:,3) = X;
            
            image(times(ix),(0:size(M,1)-1)/4,M)
            colormap((1-gray).^2)
            hold on
            
            for i=traceIdx
                times = (0:size(aodTraces(i).trace,1)-1)/aodTraces(i).fs;
                ix = times>=r(1) & times<=r(2);
                y = i-traceIdx(1);
                plot(times(ix), -aodTraces(i).trace(ix)/step+y,'Color',[0 0.2 0],'LineWidth',0.25)
                %plot(times(ix), y,':','Color',[0 0 0],'LineWidth',0.25)
            end
            hold off
            set(gca,'YTick',traceIdx-traceIdx(1),'YTickLabel',traceIdx)
            set(gca,'XTick',[10 30 50])
            ylim([-0.5 length(traceIdx)])
            ylabel 'cell number'
            xlabel 'time (s)'
            
            fig.cleanup
            fig.save('~/cov/figures/src/aod-traces.eps')
        end
        
        
        function colorbar
            r = 0.1;  % range
            fig = Figure(1,'size',[15 10]);
            s = reshape((covest.doppler),[], 1, 3);
            sz = size(s,1);
            s = permute(s,[2 1 3]);
            image(s)
            set(gca,'YTick',[],'XTick',[1 (sz+1)/2 sz], 'XTickLabel', [-r 0 r])
            set(gca,'XAxisLocation','bottom')
            set(gca,'Position',[0.1 0.6 0.8 0.25])
            fig.cleanup
            fig.save('~/cov/figures/src/colorbar-10-horiz.eps')
        end
        
        
        function aodCov
            key = fetch(pop.AodBinnedTraces2 & 'mod(aod_scan_start_time,10000)=4328' & 'bin_opt=0');
            fig = Figure(1,'size',[35 32]);
            
            C0 = fetch1(pop.AodCovEstimate2 & key & 'cov_estim_num=0','cov_matrix');  % noise correlation rather than plain correlation
            C = C0;
            p = size(C,1);
            C = corrcov(C);
            imagesc(C,.1*[-1 1]);
            colormap(covest.doppler)
            axis image
            set(gcf,'PaperPosition', [0 0 4 2.8], 'PaperSize', [4 2.8])
            set(gca,'YTick',[1 size(C,1)],'XTick',[])
            set(gca, 'Position', [0.20 0.04 0.795 0.89]);
            fig.cleanup
            box on
            fig.save('~/cov/figures/src/aod-corr.eps')
            
            close all
            delete(fig)
            fig = Figure(1,'size',[35 30]);
            [x,y] = meshgrid(1:p,1:p);
            [a,b] = hist(C(x<y),-0.6:0.002:0.2);
            area(b,a,'LineWidth',0.5,'EdgeColor',[0 0 0.4],'FaceColor',[0.7 0.7 0.8]);
            hold on
            plot([0 0],[0 1]*max(a)*1.1,'k')
            plot([1 1]*mean(C(x<y)),[0 1]*max(a)*1.1,'r-','LineWidth',1)
            xlim([-0.06 0.18])
            ylim([0 1.1*max(a)])
            xlabel correlations
            set(gca,'YTick',[])
            set(gca,'XTick',[-1:1]*.05)
            set(gca,'Position',[0.1 .3 1.0 0.7])
            
            set(gca,'YColor',[1 1 1]*.999,'TickDir','out')
            set(gcf,'PaperPosition', [0 0 2 2], 'PaperSize', [2 2])
            fig.cleanup
            fig.save('~/cov/figures/src/corr-hist.eps')
        end
        
        
        function recon
            key = fetch(pop.AodBinnedTraces2 & 'mod(aod_scan_start_time,10000)=4328' & 'bin_opt=0');
            C0 = fetch1(pop.AodCovEstimate2 & key & 'cov_estim_num=0','cov_matrix');  % noise correlation rather than plain correlation
            cove = fetch(pop.AodCovEstimate2 & key & 'cov_estim_num=14','*');  % noise correlation rather than plain correlation
            p = size(C0,1);
            CC0 = corrcov(C0);
            
            % correlation matrix: sample / regularized
            comboPlot(CC0,corrcov(cove.cov_matrix),...
                '~/cov/figures/src/combo-corrs.eps')
            densityPlot(corrcov(C0),corrcov(cove.cov_matrix),...
                '~/cov/figures/src/scatter-corrs')
            
            % partial correlation matrix: sample / regularized
            iC0 = -corrcov(inv(C0));
            iC0 = (~eye(p)).*iC0;
            iC = inv(cove.cov_matrix);
            %ds = diag(sqrt(diag(iC)));
            ds = diag(sqrt(diag(cove.sparse)));
            iC = -(~eye(p)).*(ds\iC/ds);
            
            % partial correlation matrix: sample / regularized
            comboPlot(iC0,iC,'~/cov/figures/src/combo-pcorrs.eps')
            densityPlot(iC0,iC,'~/cov/figures/src/scatter-pcorrs')
            % replace interaction matrix with thresholded correlations
            [ii,jj] = meshgrid(1:p,1:p);
            threshold = quantile(abs(CC0(ii>jj)),covest.sparsity(cove.sparse));
            
            densityPlot(CC0,-ds\cove.sparse/ds,...
                '~/cov/figures/src/scatter-corrs-vs-sparse',[0 0.2],threshold)
            
            comboPlot(-ds\cove.sparse/ds,ds\cove.lowrank*cove.lowrank'/ds,...
                '~/cov/figures/src/combo-SL.eps')
            
            
            function densityPlot(C1,C2,filename,axisAlphas,threshold)
                if nargin<4
                    axisAlphas=[.2 .2];
                end
                rng=.1;
                clf
                pp = size(C1,1);
                [i,j] = ndgrid(1:pp,1:pp);
                C1 = double(C1(i>j));
                C2 = double(C2(i>j));
                r = [-.7 1.9]*rng;
                i = linspace(r(1),r(2),200);
                [~,ix] = min(abs(i));
                i = i - i(ix);
                r = [min(i) max(i)];
                C1 = interp1(i,1:length(i),max(r(1),min(r(2),C1)),'nearest');
                C2 = interp1(i,1:length(i),max(r(1),min(r(2),C2)),'nearest');
                I = accumarray([C2 C1],1,length(i)*[1 1]);
                cmap=[[1 1 1];jet(100)];
                I = reshape(cmap(min(I+1,size(cmap,1)),:),[size(I) 3]);
                if nargin>=5
                    I(:,:,:) = bsxfun(@times,...
                        I(:,:,:),.8+.2*double(i'<1000)*double((abs(i)>=threshold)));
                end
                image(i,i,I)
                set(gca,'XTick',rng,'YTick',rng)
                set(refline(1,0),'Color','r','LineWidth',.5)
                set(gca,'Position',[0.02 0.02 0.96 0.96])
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
                colormap(covest.doppler)
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
        
        
        function network(scanNum)
            showInteractions = true;
            alpha = 0.05;
            if ~nargin
                scanNum = 4328;
            end
            
            radius = 200;
            paperSize = [12 11];
            xticks = (-1:1)*50;
            zticks = [200 250];
            zm = .6;
            panx = 1.3;
            pany = -.8;
            alphaMultiplier = 2.0;
            lineWidth = .5;
            yoffset = -15;
            xoffset = 10;
            zoffset = 25;
            
            filename = sprintf('~/cov/figures/src/network-%04d',scanNum);
            doCorrs = false;
            if true
                filename = sprintf('~/cov/figures/src/network-%04d-restricted',scanNum);
                radius = 45;  % limiting radius.  Make large (>300 to plot all)
                paperSize = [6.5 5.0];
                xticks = (-2:2)*25;
                zticks = [200 240];
                zm = .73;
                panx = 1.1;
                pany = -2.0;
                alphaMultiplier = 2;
                lineWidth = 1;
            end
            
            close all
            
            key = fetch(pop.AodCovEstimate2 & ...
                sprintf('mod(aod_scan_start_time,10000)=%d',scanNum) & ...
                'bin_opt=0' & 'cov_estim_num=14');
            
            [S,cellnums] = fetch1(pop.AodCovEstimate2*pop.AodBinnedTraces2 & key,...
                'sparse','cellnums');
            S = corrcov(S);
            C0 = fetch1(pop.AodCovEstimate2 & setfield(key,'cov_estim_num',0), 'cov_matrix');
            C0 = corrcov(C0);
            
            if doCorrs
                disp 'remove this  block!!!'
                p = size(S,1);
                % replace interaction matrix with thresholded correlations
                [ii,jj] = meshgrid(1:p,1:p);
                threshold = quantile(abs(C0(ii>jj)),covest.sparsity(S));
                S = -C0.*(abs(C0)>threshold);
                filename = sprintf('~/cov/figures/src/network-%04d-thresh-corrs',scanNum);
            end
            [ori,pval,oriCellnums] = fetch1(pop.AodVonMises & key, ...
                'von_pref', 'von_p_value', 'cellnums');
            
            [x,y,z,xyzcellnums] = fetchn(aod.Traces*aod.UniqueCell & key, ...
                'x','y','z','cell_num');
            select = (x-xoffset).^2+(y-yoffset).^2+(z-zoffset).^2 < radius.^2 & ismember(xyzcellnums,cellnums);
            assert(all(oriCellnums(:)==xyzcellnums(:)))
            x = x(select);
            y = y(select);
            z = z(select);
            xyzcellnums = xyzcellnums(select);
            ori = ori(select);
            pval = pval(select);
            clear oriCellnums
            z = z+200;
            ix = ismember(cellnums,xyzcellnums);
            S = S(ix,ix);
            
            % plot balls
            hue = mod(ori(:)/pi,1);
            sat = pval(:)<alpha;
            val = 1-.2*(pval(:)>=alpha);
            color = hsv2rgb([hue sat val]);
            
            ne7.vis.scatter3sph(x,y,z,'siz',3,'col',color)
            light('Position',[0.5 0.5 1],'Style','infinit','Color',[1 1 1])
            lighting gouraud
            axis vis3d
            axis equal
            set(gca,'ZDir','reverse')
            camproj perspective
            grid on
            
            % add interactions
            if showInteractions
                [i,j] = ndgrid(1:size(S,1),1:size(S,2));
                % positive interactions
                ix = find(j(:)>i(:) & S(:)<0);
                for ix = ix'
                    lineAlpha = min(1,(alphaMultiplier*1*abs(S(ix)))).^.9;
                    patch([x(i(ix)) x(j(ix))],[y(i(ix)) y(j(ix))], [z(i(ix)) z(j(ix))], [0 0], ...
                        'EdgeColor',[0 1 0]/2,'EdgeAlpha',lineAlpha,'LineWidth',lineWidth,'FaceColor','none')
                end
                % negative interactions
                ix = find(j(:)>i(:) & S(:)>0);
                % ix = [];
                for ix = ix'
                    lineAlpha = min(1,(alphaMultiplier*1.5*abs(S(ix)))).^.9;
                    patch([x(i(ix)) x(j(ix))],[y(i(ix)) y(j(ix))], [z(i(ix)) z(j(ix))], [0 0], ...
                        'EdgeColor',[1 0 1]/2,'EdgeAlpha',lineAlpha,'LineWidth',lineWidth,'FaceColor','none')
                end
            end
            view(15, 40)
            zoom(zm)
            camPos = get(gca,'CameraPosition');
            set(gca,'CameraPosition',camPos*0.5)
            set(gca,'ZTick',zticks)
            set(gca,'XTick',xticks,'YTick',xticks)
            campan(panx,pany)
            set(gca,'fontsize',8,'linewidth',0.25,'TickLength',get(gca,'TickLength')*0.75)
            
            set(gcf,'PaperUnits','centimeters','PaperSize',paperSize,'PaperPosition',[0 0 paperSize])
            print('-dpdf','-r800',filename)
        end
        
        
        function networkMovie(scanNum)
            if ~nargin
                scanNum = 4328;
            end
            alphaMultiplier = 1.0;
            lineWidth = 1;
            alpha = 0.05;
            key = fetch(pop.AodCovEstimate2 & ...
                sprintf('mod(aod_scan_start_time,10000)=%d',scanNum) & ...
                'bin_opt=0' & 'cov_estim_num=14');
            
            close all
            
            disp(key)
            writer = VideoWriter(sprintf('scan%05d-neg',scanNum),'MPEG-4');
            writer.Quality = 100;
            writer.FrameRate = 30;
            writer.open
            
            S = fetch1(pop.AodCovEstimate2 & key,'sparse');
            cellnums = fetch1(pop.AodBinnedTraces2 & key, 'cellnums');
            fprintf('%d neurons \n', length(cellnums));
            [ori,pval] = fetch1(pop.AodVonMises & key, 'von_pref', 'von_p_value');
            
            [x,y,z,xyzcellnums] = fetchn(aod.Traces*aod.UniqueCell & key, ...
                'x','y','z','cell_num');
            ix = ismember(xyzcellnums,cellnums);
            assert(all(xyzcellnums(ix)'==cellnums))
            x = x(ix);
            y = y(ix);
            z = z(ix);
            ori =  ori(ix);
            pval = pval(ix);
            
            % plot balls
            hue = mod(ori(:)/pi,1);
            sat = pval(:)<alpha;
            val = 1-0.7*(pval(:)>=alpha);
            color = hsv2rgb([hue sat val]);
            
            f = figure('Visible','off','Position',[0 0 800 600],'Color',[0 0 0.2]);
            light('Position',[0.5 0.5 1],'Style','infinit','Color',[0.5 0.5 0.6])
            lighting gouraud
            
            ne7.vis.scatter3sph(x,y,z,'siz',3,'col',color)
            axis vis3d
            axis equal
            set(gca,'Color','none')
            set(gca ,'XColor',[0.2 0.3 0.3])
            set(gca ,'YColor',[0.3 0.2 0.3])
            set(gca ,'ZColor',[0.3 0.3 0.2])
            camproj perspective
            grid on
            xlabel 'x (\mum)'
            ylabel 'y (\mum)'
            zlabel 'z (\mum)'
            
            % add interactions
            
            
            [i,j] = ndgrid(1:size(S,1),1:size(S,2));
            % positive interactions
            ix = find(j(:)>i(:) & S(:)<0);
            ix = [];
            for ix = ix'
                lineAlpha = min(1,(alphaMultiplier*3*abs(S(ix)))).^.8;
                patch([x(i(ix)) x(j(ix))],[y(i(ix)) y(j(ix))], [z(i(ix)) z(j(ix))], [0 0], ...
                    'EdgeColor',[0 1 0],'EdgeAlpha',lineAlpha,'LineWidth',lineWidth,'FaceColor','none')
            end
            % negative interactions
            ix = find(j(:)>i(:) & S(:)>0);
            for ix = ix'
                lineAlpha = min(1,(alphaMultiplier*5*abs(S(ix)))).^.8;
                patch([x(i(ix)) x(j(ix))],[y(i(ix)) y(j(ix))], [z(i(ix)) z(j(ix))], [0 0], ...
                    'EdgeColor',[1 0 1],'EdgeAlpha',lineAlpha,'LineWidth',lineWidth,'FaceColor','none')
            end
            
            last = -300;
            for i=-180:0.5:180
                if last+15<i
                    fprintf .
                    last = i;
                end
                view(45+cos(i/180*pi)*90,30)
                set(gca,'CameraPosition',get(gca,'CameraPosition').*(0.7+[0.15 0.15 0.3].*cos((i/90+[-0.3 0.0 0.3])*pi)))
                writer.writeVideo(imcapture(f))
            end
            fprintf \n
            writer.close
            close(f)
        end
        
        
        
        function Segment(varargin)
            
            for key = fetch(pop.Segment & varargin)'
                [g,r] = fetch1(tp.FineAlign & key, 'fine_green_img', 'fine_red_img');
                imshowpair(g,r)
                axis on
                grid on
                masks = fetchn(pop.Trace & key, 'mask_pixels');
                
                if ~isempty(masks)
                    bw = false(size(g));
                    for m = masks'
                        bw(m{:}) = true;
                    end
                    hold on
                    b = bwboundaries(bw,4);
                    for i=1:length(b)
                        plot(b{i}(:,2),b{i}(:,1),'r')
                    end
                    hold off
                end
            end
        end
        
        
        
        
        function frac = OriChange(varargin)
            count = 0;
            ids = [];
            frac = [];
            for key = fetch(pop.GraphicalModel*pop.ModelOpt & varargin &  'pairwise')'
                ids = union(ids,key.animal_id);
                s1 = fetch(pop.GraphicalModel & key, '*');
                key.model_opt = 1;
                s0 = fetch(pop.GraphicalModel & key, '*');
                s1.stim_drives = s1.stim_drives(1:size(s0.stim_drives,1),:);
                ix = any(s1.stim_drives,2) & any(s0.stim_drives,2);
                count = count + ~isempty(ix);
                ori1 = atan2(s1.stim_drives(ix,1),s1.stim_drives(ix,2))*90/pi;
                ori0 = atan2(s0.stim_drives(ix,1),s0.stim_drives(ix,2))*90/pi;
                frac(end+1,:) = [mean(any(s0.stim_drives,2)),mean(any(s1.stim_drives,2))];
            end
        end
        
        
        
        function GraphicalModel(varargin)
            
            for key = fetch(pop.GraphicalModel & varargin & 'nhidden>0')'
                clf
                [g,r] = fetch1(tp.FineAlign & key, 'fine_green_img', 'fine_red_img');
                
                subplot 231
                imshow(cat(3,r,g,0*g)/1000);
                [x,y,masks] = fetchn(pop.Trace & key, 'centroid_x','centroid_y','mask_pixels');
                outlineCells
                
                % plot connectivity with hidden unit
                subplot 233
                surf(0*g,g,'EdgeColor','none')
                colormap gray
                camproj('perspective')
                view(0,70)
                outlineCells
                s = fetch(pop.GraphicalModel & key, '*');
                if isempty(s.cell_hidden)
                    continue
                end
                axis image
                axis off
                plotConnections
                
                
                % plot hidden unit
                [sx,sy,sz]= sphere(20);
                radius = 6;
                cx = size(g,2)/2;
                cy = size(g,1)/2;
                cz = size(g,2);
                sx = sx*radius + cx;
                sy = sy*radius + cy;
                sz = sz*radius + cz;
                
                surf(sx, sy, sz,'FaceColor',[0.1 0.3 0.8],'EdgeColor','none');
                for i=1:length(s.cell_hidden(:,1))
                    if s.cell_hidden(i,1)
                        plot3([cx x(i)],[cy y(i)],[cz 0], ':', 'Color', [0.1 0.5 0.1]);
                    end
                end
                set(gca,'YDir','reverse')
                subplot 236
                imshow(g/1000)
                outlineCells
                plotOri
                
                savedKey = key;
                
                % plot raw tuning
                subplot 234
                key.model_opt = 2;
                s = fetch(pop.GraphicalModel & key, '*');
                imshow(g/1000)
                outlineCells
                plotOri
                
                % plot full-model tuning
                key.model_opt = 3;
                subplot 235
                s = fetch(pop.GraphicalModel & key, '*');
                imshow(g/1000)
                outlineCells
                plotOri
                
                % plot full-model connectivity
                subplot 232
                imshow(g/1000);
                outlineCells
                plotConnections
                
                
                set(gcf,'PaperSize',[18 9],'PaperPosition',[0 0 18 9])
                fname = sprintf('figures/pairwise_%5d_%03d_opt%u',key.animal_id, key.scan_idx, savedKey.model_opt);
                print('-dpdf','-r600', fname)
            end
            
            
            function plotOri
                for ii=1:length(masks)
                    if any(s.stim_drives(ii,:))
                        ori = atan2(s.stim_drives(ii,2),s.stim_drives(ii,1))/2;
                        hue = mod(ori/pi+0.5,1);
                        col = hsv2rgb([hue,1,1]);
                        rr = sqrt(length(masks{ii}))/2;
                        hold on
                        plot(x(ii)-sin(ori)*rr*[-1 1],y(ii)+cos(ori)*rr*[-1 1],'color',col,'LineWidth',3)
                    end
                end
            end
            
            
            function plotConnections
                K = triToSquare(s.cell_cell);
                hold on
                for ii=2:size(K,1)
                    for jj=1:ii-1
                        if K(ii,jj) < 0
                            plot(x([ii jj]),y([ii jj]),'g');
                        elseif K(ii,jj) > 0
                            plot(x([ii jj]),y([ii jj]), 'r');
                        end
                    end
                end
            end
            
            
            function outlineCells
                bw = false(size(g));
                for m = masks'
                    bw(m{:}) = true;
                end
                hold on
                b = bwboundaries(bw,4);
                for iMask=1:length(b)
                    plot(b{iMask}(:,2),b{iMask}(:,1),'Color',[0.4 0.4 1.0])
                end
                hold off
            end
        end
    end
end