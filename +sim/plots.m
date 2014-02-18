classdef plots
    methods(Static)
        function figure
            restriction = 'true_cov_id in (8,9,10,11)';
            letters = {
                'independent'   'A'
                'multifactor'    'B'
                'sparse'         'C'
                'sparse+lowrank' 'D'
                };
            
            iBasis = 1;
            for key = fetch(sim.TrueCov & restriction,'covclass','p')'
                clf
                letter = letters{strcmp(key.covclass,letters(:,1)),2};
                d = 16;  % resize factor 
                
                %%%%% PANEL 2: True covariance matrix
                Sigma = fetch1(sim.TrueCov & key, 'covmatrix');
                Sigma = corrcov(Sigma);
                cmap = covest.doppler(256);
                r = 0.25;
                Sigma = max(1,min(256,ceil(Sigma/(2*r)*256+128))); 
                Sigma = reshape(cmap(Sigma,:),[size(Sigma) 3]);
                Sigma = imresize(Sigma,d,'nearest');
                imwrite(Sigma, sprintf('~/cov/figures/src/sim-truth-%c.png',letter));
                
               
                %%%%% ROW 3: Sample covariance matrix
                sampleSize = 1000;
                C = fetch1(sim.Sample & setfield(key,'sample_size',sampleSize) & 'seed=1', 'sample_covmatrix'); %#ok<SFLD>
                C = corrcov(C);
                C = max(1,min(256,ceil(C/(2*r)*256+128))); 
                C = reshape(cmap(C,:),[size(C) 3]);
                C = imresize(C,d,'nearest');
                imwrite(C, sprintf('~/cov/figures/src/sim-sample%d-%c.png',sampleSize,letter))
                
                %%%%%%% ROW 4: Matching reconstruction
                iBasis = iBasis + 1;
                estimators = {'sample','shrink','multifactor','sparse','sparse+lowrank'};
                
                C = fetch1(sim.CovEstimate & setfield(key,'sample_size',sampleSize) ...
                    & (sim.CovEstimator & struct('estimator_name',estimators{iBasis})) & 'seed=1', 'cov_matrix'); %#ok<SFLD>
                C = corrcov(C);
                C = max(1,min(256,ceil(C/(2*r)*256+128))); 
                C = reshape(cmap(C,:),[size(C) 3]);
                C = imresize(C,d,'nearest');
                imwrite(C, sprintf('~/cov/figures/src/sim-matching%d-%c.png',sampleSize,letter))
                
                %%%% ROW 5: Excess loss
                fig = Figure(1,'size',[35 55]);
                styles = {'.', 's', 'v', 'o', '*'};
                lineStyles = {'-','-','-','-','-'};
                colors = {'k', [.0 .5 .0], [.8 0 0], [0 .4 0.9], [.4 0 .9]};
                
                for i=length(estimators):-1:1
                    s = fetch(sim.CovEstimate & key & 'sample_size>150' & (sim.CovEstimator & struct('estimator_name',estimators{i})),...
                        'excess_loss','seed','sample_size');
                    [tab,~,sampleSizes] = dj.struct.tabulate(s, 'excess_loss', 'seed', 'sample_size');
                    x = (1:length(sampleSizes)) + (i*0.04-0.08);
                    y = mean(tab);
                    e = std(tab);%/sqrt(size(tab,1));
                    plot(x,y,lineStyles{i},'Color',colors{i},'LineWidth',.5)
                    hold all
                    errorbar(x,y,e,styles{i},'Color',colors{i},'MarkerSize',4,'LineWidth',.5)
                end
                hold off
                set(gca,'XTick',1:length(sampleSizes),'XTickLabel',{'.25' '.5' '1' '2' '4'})
                set(gca,'YTick',[0 .03],'YTickLabel',{'0', '.03'})
                set(gca,'Position',[0.20 0.18 0.78 0.75])
                set(gca,'YAxisLocation','left')
                ylim([0 .032])
                xlim([0.7 5.3])
                xlabel('sample size (\times1000)')
                fig.cleanup
                fig.save(sprintf('~/cov/figures/src/excess-loss-%c.eps',letter))

                %%%%% ROW 6: relative validation loss
                fig = Figure(1,'size',[35 55]);
                s = fetch(sim.CovEstimateFold & key & 'sample_size>150' & (sim.CovEstimator & struct('estimator_name',estimators{iBasis})),...
                    'fold','cross_loss','seed','sample_size');
                basis = dj.struct.tabulate(s, 'cross_loss','fold','seed', 'sample_size');

                for i=length(estimators):-1:1
                    s = fetch(sim.CovEstimateFold & key & 'sample_size>150' & (sim.CovEstimator & struct('estimator_name',estimators{i})),...
                        'fold','cross_loss','seed','sample_size');
                    [tab,~,~,sampleSizes] = dj.struct.tabulate(s, 'cross_loss','fold','seed', 'sample_size');
                    tab = squeeze(mean(tab-basis));
                    x = (1:length(sampleSizes)) + (i*0.04-0.08);
                    y = mean(tab);
                    e = std(tab);%/sqrt(size(tab,1));
                    plot(x,y,lineStyles{i},'Color',colors{i},'LineWidth',.5)
                    hold all
                    errorbar(x,y,e,styles{i},'Color',colors{i},'MarkerSize',4,'LineWidth',.5)
                end
                hold off
                xlabel 'sample size'
                set(gca,'XTick',1:length(sampleSizes),'XTickLabel',{'.25' '.5' '1' '2' '4'})
                set(gca,'YTick',-0.01:0.01:0.01,'YTickLabel',{'-.01','0','.01'})
                set(gca,'Position',[0.20 0.18 0.78 0.75])
                set(gca,'YAxisLocation','left')
                ylim([-.011 0.011])
                xlim([0.7 5.3])
                xlabel('sample size (\times1000)')

                fig.cleanup
                fig.save(sprintf('~/cov/figures/src/validation-loss-%c.eps',letter))
            end
            
            % plot legend
            fig = Figure(1,'size',[30 40]);
            for i=1:length(estimators)
                x = [-1 0 1];
                y = -[i i i]*0.2;
                e = [.03 .01 0];
                plot(x,y,lineStyles{i},'Color',colors{i},'LineWidth',.5)
                hold on
                errorbar(x,y,e,styles{i},'Color',colors{i},'MarkerSize',4,'LineWidth',.5)
            end
            hold off
            axis off
            ylim([-1.2 0])
            xlim([-1 1]*1.05)
            set(gca,'Position',[.3 -.1 .3 1.2])
            fig.cleanup
            fig.save('~/cov/figures/src/sim-legend.eps')
        end

    end
end