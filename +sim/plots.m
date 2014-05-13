classdef plots
    properties(Constant)
        figPath = '~/cov/figures/src'
    end
    methods(Static)
        function figure1
            r = 0.25;
            cmap = uint8(cove.doppler(256)*255);
            key.truth_seed = 1;
            key.sample_seed = 1;
            key.sample_size = 1000;
            
            % Row 2
            [C,truthKeys] = fetchn(sim.Truth & 'truth_seed=1','true_cov', 'ORDER BY truth');
            for i=1:length(C)
                c = max(0,min(1,(corrcov(C{i})/r+1)/2))*255+1;
                img = reshape(cmap(round(c),:),size(c,1),size(c,1),3);
                img = imresize(img,12,'nearest'); % upsample for better rendering in PDF
                imwrite(img,fullfile(sim.plots.figPath,sprintf('Fig1-2-%c.png',64+i)))
            end
            
            % Row 3
            C = arrayfun(@(k) ...
                fetch1(sim.CovMatrix & key & k & 'method=0' & 'nfolds=1', 'cov_matrix'), ...
                truthKeys,'uni',false);
            for i=1:length(C)
                c = max(0,min(1,(corrcov(C{i})/r+1)/2))*255+1;
                img = reshape(cmap(round(c),:),size(c,1),size(c,1),3);
                img = imresize(img,12,'nearest'); % upsample
                imwrite(img,fullfile(sim.plots.figPath,sprintf('Fig1-3-%c.png',64+i)))
            end
            
            % Row 4
            methods = [10 30 80 90];
            C = arrayfun(@(k,m) ...
                fetch1(sim.CovMatrix & key & k & struct('method',m,'nfolds',1), 'cov_matrix'), ...
                truthKeys,methods','uni',false);
            for i=1:length(C)
                c = max(0,min(1,(corrcov(C{i})/r+1)/2))*255+1;
                img = reshape(cmap(round(c),:),size(c,1),size(c,1),3);
                img = imresize(img,12,'nearest'); % upsample
                imwrite(img,fullfile(sim.plots.figPath,sprintf('Fig1-4-%c.png',64+i)))
            end
            clear all
            
            % Row 5
            s = fetch(sim.CovMatrix*sim.TrueLoss*pro(sim.Truth,'truth_type') & 'nfolds=1', 'truth_type', 'true_loss', 'ORDER BY truth, method, sample_size');
            [loss,truth,method,sampleSize] = dj.struct.tabulate(s,'true_loss','truth_type','method','sample_size');
            xticks = sampleSize/1000;
            yticks = [0 0.04];
            nSamples = size(s,4);
            offset = exp(linspace(-.1,.1,length(method)));
            for i=1:length(truth)
                fig = Figure(1,'size',[35 55]);
                m = squeeze(mean(loss(i,:,:,:),4))';
                e = squeeze(std(loss(i,:,:,:),[],4))'/sqrt(nSamples);
                
                x = exp(bsxfun(@plus, log(sampleSize'), 0.05*linspace(-1,1,length(method))));
                h = plot(x, m);
                hold on
                errorbar(x,m,e);
                set(gca,'XTick',sampleSize,'XTickLabel',nozero(xticks),...
                    'YTick',yticks,'YTickLabel',nozero(yticks), ...
                    'XScale','log','Position',[.2 .1 .8 .9]);
                ylim([0 .05])
                xlim([200 5000])
                fig.cleanup
                fig.save(fullfile(sim.plots.figPath, sprintf('Fig1-5-%c.eps',64+i)))
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