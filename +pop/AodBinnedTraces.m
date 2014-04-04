%{
pop.AodBinnedTraces (computed) # aod traces binned for correlation computations
-> pop.StableScans
-> aod.UniqueCells
-> aod.TracePreprocessSet
-> pop.BinOpt
---
nneurons                    : smallint                      # the number of neurons included in the analysis
cellnums                    : longblob                      # list of participating cell numbers
ndirs=null                  : smallint                      # number of directions
dirs=null                   : longblob                      # direction of grating movement
psths=null                  : longblob                      # peristimulus histograms for each cell and direction
binned_traces               : longblob                      # array of binned traces
binsize_ms                  : float                         # true bin size
epoch_seconds               : float                         # the duration of the analyzed episode
%}

classdef AodBinnedTraces < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = (aod.TracePreprocessSet*aod.UniqueCells*pop.StableScans*pop.BinOpt) ...
            & acq.AodStimulationLink & 'preprocess_method_num=5' & 'bin_opt=0'
    end
    
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            disp 'loading traces...'
            opt = fetch(pop.BinOpt&key,'*');
            [X,keys] = fetchn(aod.TracePreprocess & aod.UniqueCell & key, 'trace');
            
            X = [X{:}];
            v = var(X);
            thresh = opt.activity_thresh*quantile(v,0.9);
            ix = var(X)>thresh;
            fprintf('found %d active cells of %d total\n',sum(ix),size(X,2))
            X = X(:,ix);
            keys = keys(ix);
            key.cellnums = [keys.cell_num];
            key.nneurons = length(key.cellnums);
            
            % subtract mean
            M = mean(X);
            X = bsxfun(@minus, X, M);
            
            % subtract stimulus
            if opt.psth_ms>0
                disp 'subtracting stimulus response...'
                [S,psth,oris] = pop.AodBinnedTraces2.computeStim(key,X,0,opt.psth_ms);
                key.ndirs = length(oris);
                key.dirs = oris;
                key.psths = psth;
                
                % limit analysis to stimulus period
                ix = find(any(S,2));
                ix = ix(1):ix(end);
                X = X(ix,:);
                X = bsxfun(@minus,X,mean(X));
                S = S(ix,:);
                S = bsxfun(@minus, S, mean(S));  % center stimulus response on zero
                X = X - S;  % subtract stimulus
            end
            
            disp binnng..
            fps = 1000/median(diff(getTimes(aod.TracePreprocess & key)));
            key.epoch_seconds = size(X,1)/fps;
            framesPerBin = max(1,round(opt.bin_ms/1000*fps));
            
            if framesPerBin
                % rebin
                X = X(1:floor(end/framesPerBin)*framesPerBin,:);
                Y = 0;
                for offset=1:framesPerBin
                    Y = Y + X(offset-1+(1:framesPerBin:end),:)/framesPerBin;
                end
                X = Y;
            end
            key.binsize_ms = key.epoch_seconds/size(X,1)*1000;
            key.binned_traces = single(X);
            self.insert(key)
        end
    end
    
    
    
    methods (Static)
        
        function [S,psth,oris] = computeStim(key,X,lag,duration)
            % compute the average stimulus trace from traces
            
            times = getTimes(aod.TracePreprocess & key);
            
            % limit the trial group to 2-stimulus trials
            trialGroup = stimulation.StimTrialGroup*acq.AodStimulationLink;
            trialGroup = pro(trialGroup, stimulation.StimConditions, 'count(*)->ncond') & key & 'ncond<8';
            assert(trialGroup.count == 1);
            trials = fetch(stimulation.StimTrials*trialGroup);
            conditions = fetch(stimulation.StimConditions*trialGroup,'*');
            
            oris = unique(arrayfun(@(x) x.orientation, [conditions.condition_info]));
            framesInTrial = floor(duration/median(diff(times)));
            
            % Extract segment of trials for each stimulus
            presentations = [];
            for i = 1:length(trials)
                trial_info = fetch1(stimulation.StimTrials & trials(i), 'trial_params');
                event = fetch(stimulation.StimTrialEvents(trials(i), 'event_type="showSubStimulus"'),'*');
                onset = sort([event.event_time]);
                for j = 1:length(event)
                    cond = trial_info.conditions(j);
                    ori = conditions(cond).condition_info.orientation;
                    condIdx = find(ori == oris);
                    idx = find(times >= (onset(j) + lag),1,'first');
                    if ~isempty(idx) && idx<=size(X,1)-framesInTrial
                        s.idx = idx;
                        s.ori = ori;
                        s.condIdx = condIdx;
                        s.traces = X(idx+(0:framesInTrial-1),:);
                        presentations = [presentations; s]; %#ok<AGROW>
                    end
                end
            end
            
            % compute PSTHs
            for iCond = 1:length(oris)
                t = cat(3, presentations([presentations.condIdx]==iCond).traces);
                psth{iCond} = mean(t,3); %#ok<AGROW>
            end
            
            % compute response traces
            S = zeros(size(X));
            for p = presentations'
                S(p.idx+(0:framesInTrial-1),:) = ...
                    S(p.idx+(0:framesInTrial-1),:) + psth{p.condIdx};
            end
        end
    end
end
