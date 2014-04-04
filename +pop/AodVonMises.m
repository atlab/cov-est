%{
pop.AodVonMises (computed) # von Mises tuning for all cells
-> aod.TracePreprocessSet
-> aod.UniqueCells
-----
cellnums     : blob      # cells used
responses    : longblob  # responses of cells to orientations
nshuffles    : smallint  # number of shuffles
von_p_value  : blob      # von mises shuffle p-values
von_r2       : blob      # R-squared of response
von_pref     : blob      # von mises preferred direction
von_base     : blob      # von mises bases
von_amp1     : blob      # von mises preferred amplitude
von_amp2     : blob      # von mises anti-preferred amplitude
von_sharp    : blob      # von mises sharpness
%}

classdef AodVonMises < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = (aod.TracePreprocessSet*aod.UniqueCells*pop.StableScans) ...
            & acq.AodStimulationLink & 'preprocess_method_num=5'
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            disp 'loading traces...'
            [X,cellnums] = fetchn(aod.TracePreprocess & aod.UniqueCell & key, ...
                'trace', 'cell_num', 'ORDER BY cell_num');
            X = [X{:}];
            times = getTimes(aod.TracePreprocess & key);
            assert(size(X,1)==length(times))
            
            % limit the trial group to 2-stimulus trials
            trialGroup = stimulation.StimTrialGroup*acq.AodStimulationLink;
            trialGroup = pro(trialGroup, stimulation.StimConditions, ...
                'count(*)->ncond') & key & 'ncond>=8';
            assert(trialGroup.count == 1)
            trials = fetch(stimulation.StimTrials*trialGroup);
            conditions = fetch(stimulation.StimConditions*trialGroup,'*');
            
            lag = 20;  % delay (ms)
            duration = 450; % integration duration (ms)
            dt = median(diff(times));
            framesPerTrial = ceil(duration/dt);
            
            disp 'extracting events...'
            presentations = [];
            cellnums = num2cell(cellnums);
            for i = 1:length(trials)
                trial_info = fetch1(stimulation.StimTrials & trials(i), 'trial_params');
                event = fetch(stimulation.StimTrialEvents(trials(i), 'event_type="showSubStimulus"'),'*');
                onset = sort([event.event_time]);
                
                for j = 1:length(event)
                    cond = trial_info.conditions(j);
                    ori = conditions(cond).condition_info.orientation;
                    
                    idx = find(times >= (onset(j) + lag),1,'first');
                    if ~isempty(idx) && idx<=size(X,1)-framesPerTrial
                        responses = num2cell(mean(X(idx+(0:framesPerTrial-1),:)));
                        presentations = [presentations;
                            struct('ori',ori,'cellnum',cellnums,'response',responses')]; %#ok<AGROW>
                    end
                end
            end
            disp 'tabulating...'
            [responses,cellnums,oris] = dj.struct.tabulate(presentations, 'response', 'cellnum', 'ori');
            disp 'computing von Mises tuning...'
            
            nShuffles = 10000;
            [von, r2, p] = ne7.rf.VonMises2.computeSignificance(responses, nShuffles);
            key.cellnums = int16(cellnums);
            key.responses = single(responses);
            key.nshuffles = nShuffles;
            key.von_r2 = single(r2);
            key.von_base = single(von.w(:,1));
            key.von_amp1 = single(von.w(:,2));
            key.von_amp2 = single(von.w(:,3));
            key.von_sharp= single(von.w(:,4));
            key.von_pref = single(von.w(:,5));
            key.von_p_value = single(p);
            
            self.insert(key)
        end
    end
end
