% plots and figures for the covariance estimation paper


classdef plots
    
    methods(Static)
        function compareCorrMethods
            % compare correlation coefficients with variance estimated in
            % each bin or globally
            
            for key = fetch(covest.ActiveCells & 'preprocess_method_num=5' & 'high_repeats')'
                [X,evokedBins,ix] = fetch1(covest.ActiveCells*covest.Traces & key, ...
                    'trace_segments','evoked_bins', 'selection');
                evokedBins = min(evokedBins,size(X,1));
                
                % compare with original binned traces
                ix = ismember(...
                    fetchn(aod.TracePreprocess & aod.UniqueCell & key, 'cell_num'), ...
                    fetch1(pop.AodBinnedTraces & rmfield(key,'preprocess_method_num'), 'cellnums'));
                X = X(1:min(end,evokedBins),:,:,ix);
                %X = X(:,:,:,ix);
                
                Y = fetch1(pop.AodBinnedTraces & rmfield(key,'preprocess_method_num'), 'binned_traces');
                
                [nBins,nConds,nTrials,nCells] = size(X);
                
                % subtract mean response
                if size(X,1)<=evokedBins
                    M = nanmean(X,3);
                else
                    M1 = nanmean(X(1:evokedBins,:,:,:),3);   % binwise mean of evoked response
                    M2 = reshape(nanmean(reshape(X(evokedBins+1:end,:,:,:),[],nCells)),1,1,1,nCells);  % common mean in intertrial periods
                    M2 = repmat(M2,size(X,1)-evokedBins,nConds);
                    M = cat(1,M1,M2);
                end
                
                X = bsxfun(@minus, X, M);
                
                % method 0: original computation
                c0 = corr(Y);
                
                % method 1: common variance
                V = reshape(nanvar(reshape(X,[],nCells)), 1,1,1,nCells);
                c1 = getCorr(X,V);
                
                % method 2: condition-specific variance
                V = reshape(nanvar(reshape(permute(X, [1 3 2 4]),[],nConds,nCells)), 1,nConds,1,nCells);
                c2 = getCorr(X,V);
                
                % method 3: bin-specific variance
                V = nanvar(X,[],3);
                c3 = getCorr(X,V);
                
                
                cc1 = c1;
                cc2 = c3;
                clear c0 c1 c2 c3
                
                
                subplot 231, imagesc(cc1,[-1 1]*.2), axis image off, title 'common variance'
                colorbar
                subplot 232, imagesc(cc2,[-1 1]*.2), axis image off, title 'binned variance'
                colorbar
                colormap(covest.lib.doppler)
                
                n = size(cc1,1);
                [i,j] = meshgrid(1:n,1:n);
                cc1 = cc1(i<j);
                cc2 = cc2(i<j);
                mm1 = mean(cc1);
                mm2 = mean(cc2);
                
                r = [-.1 .3];
                subplot 234, hist(cc1,linspace(r(1),r(2),100)), xlim(r), ylim_ =ylim; xlabel correlations, hold on, plot(mm1*[1 1],ylim_,'r'), hold off, grid on
                subplot 235, hist(cc2,linspace(r(1),r(2),100)), xlim(r), ylim(ylim_), xlabel correlations, hold on, plot(mm2*[1 1],ylim_,'r'), hold off, grid on
                
                ix = linspace(r(1)/2,r(2)/2,50);
                
                subplot (2,3,[3 6]), densityPlot(ix, ix, cc1, cc2), xlabel 'common variance', ylabel 'binned variance'
                hold on, plot(mm1,mm2,'r+','MarkerSize',40), hold off
                
                set(gcf,'PaperSize',[12,8],'PaperPosition',[0 0 12 8])
                print('-dpng','-r400', sprintf('~/dev/temp/corrComp5-%04d',mod(key.aod_scan_start_time,1e4)))
            end
            
            
            function c = getCorr(X,V)
                z = reshape(bsxfun(@rdivide,X,sqrt(V+eps)),[],size(X,4));
                c = corrcov(nancov(z));
            end
            
        end
        
    end
end


function densityPlot(x,y,c1,c2)
img = hist3([c2 c1],{x,y});
img = 255*img/max(img(:));
cmap = [1 1 1; jet(255)];
%(1-bsxfun(@power, gray(256), [.2 .3 .4])
subimage(x, y, img, cmap)
set(gca,'YDir','normal')
grid on
end