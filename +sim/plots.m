classdef plots
    properties(Constant)
        figPath = '~/cov/figures/src'
    end
    methods(Static)
        function figure1
            r = 0.25;
            cmap = uint8(cove.doppler(256)*255);
            key.truth_seed = 1;
            key.sample_size = 500;
            
            
            % Row 2
            truthTypes = 'truth_type in ("diag","factor","sparse refit","sparse+latent refit")';
            [C,truthKeys] = fetchn(sim.Truth & 'truth_seed=1' & truthTypes,'true_cov', 'ORDER BY truth');
            for i=1:length(C)
                c = max(0,min(1,(corrcov(C{i})/r+1)/2))*255+1;
                img = reshape(cmap(round(c),:),size(c,1),size(c,1),3);
                img = imresize(img,12,'nearest'); % upsample for better rendering in PDF
                imwrite(img,fullfile(sim.plots.figPath,sprintf('Fig1-2-%c.png',64+i)))
            end
            
            for m=1:2
                model.model = m;
               
                % Row 3
                C = arrayfun(@(k) ...
                    fetch1(sim.CovMatrix & key & k & 'method=0' & model & 'nfolds=1', 'corr_matrix'), ...
                    truthKeys,'uni',false);
                for i=1:length(C)
                    c = max(0,min(1,(corrcov(C{i})/r+1)/2))*255+1;
                    img = reshape(cmap(round(c),:),size(c,1),size(c,1),3);
                    img = imresize(img,12,'nearest'); % upsample
                    imwrite(img,fullfile(sim.plots.figPath,sprintf('Fig1-3-%c-%u.png',64+i,model.model)))
                end
                
                % Row 4
                methods = [10 30 80 90];
                C = arrayfun(@(k,m) ...
                    fetch1(sim.CovMatrix & key & k & struct('method',m,'nfolds',1) & model, 'corr_matrix'), ...
                    truthKeys,methods','uni',false);
                for i=1:length(C)
                    c = max(0,min(1,(corrcov(C{i})/r+1)/2))*255+1;
                    img = reshape(cmap(round(c),:),size(c,1),size(c,1),3);
                    img = imresize(img,12,'nearest'); % upsample
                    imwrite(img,fullfile(sim.plots.figPath,sprintf('Fig1-4-%c-%u.png',64+i,model.model)))
                end
                clear methods
                
                % Row 5
                styles = {'.-', 's-', 'v-', 'o-', '*-'};
                colors = {'k', [.0 .5 .0], [.8 0 0], [0 .4 0.9], [.4 0 .9]};
                s = fetch(sim.CovMatrix*sim.TrueLoss*pro(sim.Truth,'truth_type') ...
                    & 'nfolds=1' & truthTypes & model, ...
                    'truth_type', 'true_loss', 'ORDER BY truth, method, sample_size');
                [loss,truth,method,sampleSize] = ...
                    dj.struct.tabulate(s,'true_loss','truth_type','method','sample_size');
                xticks = sampleSize/1000;
                yticks = 0:0.02:0.1;
                nSamples = size(loss,4);
                offset = -0.5*linspace(-.1,.1,length(method));
                x = exp(bsxfun(@plus, log(sampleSize'), offset));
                for i=1:length(truth)
                    fig = Figure(1,'size',[35 55]);
                    m = squeeze(mean(loss(i,:,:,:),4))';
                    e = squeeze(std(loss(i,:,:,:),[],4))';%/sqrt(nSamples);
                    for iMethod=1:size(m,2)
                        plot(x(:,iMethod), m(:,iMethod), styles{iMethod}, ...
                            'Color',colors{iMethod}, 'LineWidth',.5, 'MarkerSize',4)
                        hold on
                        myerrorbar(x(:,iMethod),m(:,iMethod),e(:,iMethod), ...
                            'Marker','none', 'Color',colors{iMethod}, 'LineWidth',.5)
                    end
                    hold off
                    
                    set(gca,'XTick',sampleSize,'XTickLabel',nozero(xticks),...
                        'YTick',yticks,'YTickLabel',nozero(yticks), ...
                        'XScale','log','Position',[.2 .18 .8 .8])
                    ylim([0 .055])
                    xlim([200 5000])
                    xlabel 'sample size \times1000'
                    fig.cleanup
                    fig.save(fullfile(sim.plots.figPath, sprintf('Fig1-5-%c-%u.eps',64+i,model.model)))
                end
                
                % Row 6
                s = fetch(sim.CovMatrix*pro(sim.Truth,'truth_type') & model ...
                    & 'nfolds>1' & truthTypes, ...
                    'truth_type', 'cv_loss', 'k', 'ORDER BY truth, method, sample_size');
                [loss,truth,~,sampleSize,~] = ...
                    dj.struct.tabulate(s,'cv_loss','truth_type','method','sample_size','k');
                loss = squeeze(mean(loss,4));
                for i=1:length(truth)
                    loss(i,:,:,:) = bsxfun(@minus, loss(i,:,:,:), loss(i,i+1,:,:));
                end
                
                for i=1:length(truth)
                    fig = Figure(1,'size',[35 55]);
                    m = squeeze(mean(loss(i,:,:,:),4))';
                    e = squeeze(std(loss(i,:,:,:),[],4))'; %/sqrt(nSamples);
                    for iMethod=1:size(m,2)
                        plot(x(:,iMethod), m(:,iMethod), styles{iMethod}, ...
                            'Color',colors{iMethod}, 'LineWidth',.5, 'MarkerSize',4)
                        hold on
                        myerrorbar(x(:,iMethod),m(:,iMethod),e(:,iMethod), ...
                            'Marker','none', 'Color',colors{iMethod}, 'LineWidth',.5)
                    end
                    hold off
                    
                    set(gca,'XTick',sampleSize,'XTickLabel',nozero(xticks),...
                        'YTick',yticks,'YTickLabel',nozero(yticks), ...
                        'XScale','log','Position',[.2 .18 .8 .8])
                    ylim([-0.002 .045])
                    xlim([200 5000])
                    xlabel 'sample size \times1000'
                    fig.cleanup
                    fig.save(fullfile(sim.plots.figPath, sprintf('Fig1-6-%c-%u.eps',64+i, model.model)))
                end
            end
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


function myerrorbar(x,y,e,varargin)
assert(size(x,2)==1 && size(y,2)==1 && size(e,2)==1)
n = size(x,1);
assert(size(y,1)==n && size(e,1)==n)

for i=1:n
    line(x([i i]),y(i)+e(i)*[-1 1],varargin{:})
end

end