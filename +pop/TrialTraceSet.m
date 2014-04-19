%{
pop.TrialTraceSet (computed) # traces associated with trials, aligned on stimulus onset
-> pop.StableScans
-> aod.UniqueCells
-> aod.TracePreprocessSet
high_repeats                : tinyint                       # 1=this dataset has few stimulus conditions with many trials
---
latency_ms                  : smallint                      # (ms) presumed screen-to-V1 latency
bin_ms                      : float                         #
evoked_bins                 : tinyint                       # number of evoked bins -- the remainder are spontaneous
nneurons                    : smallint                      # the number of neurons included in the analysis
cell_xyz                    : blob                          #
cellnums                    : blob                          # cell selection
ndirs                       :  tinyint                       # number of directions
directions                  : blob                          # directions
ntrials                     : blob                          # trial in each direction
trace_segments              : longblob                      # trace segements: nBins x nDirs x nTrials x nCells
%}

classdef TrialTraceSet < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = (aod.TracePreprocessSet*aod.UniqueCells*pop.StableScans) & acq.AodStimulationLink & 'preprocess_method_num in (5,7)'
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            binDuration = 150; % (ms)
            latency = 20;  % (ms) visual latency
            
            disp 'loading traces...'
            times = getTimes(aod.TracePreprocess & key);
            [X,keys] = fetchn(aod.TracePreprocess & aod.UniqueCell & key, 'trace');
            X = [X{:}];
            
            for highRepeats = [false true]
                tuple = key;
                tuple.high_repeats = highRepeats;
                tuple.bin_ms = binDuration;
                nCells = length(keys);
                fprintf('found %d cells\n',nCells)
                tuple.nneurons = nCells;
                tuple.cellnums = [keys.cell_num];
                tuple.latency_ms = latency;
                nCells = length(tuple.cellnums);
                
                % load cell locations
                list = sprintf('%d,',tuple.cellnums);
                [x,y,z] = fetchn(aod.Traces*aod.UniqueCell & key & sprintf('cell_num in (%s)', list(1:end-1)), ...
                    'x','y','z','ORDER BY cell_num');
                tuple.cell_xyz = [x y z];
                assert(nCells == size(tuple.cell_xyz,1))
                
                % find trial groups with fewer than 8 conditions
                trialGroup = stimulation.StimTrialGroup*acq.AodStimulationLink;
                
                trialGroup = pro(trialGroup, stimulation.StimConditions, 'count(*)->ncond') & key;
                if highRepeats
                    trialGroup = trialGroup & 'ncond < 8';
                else
                    trialGroup = trialGroup & 'ncond >= 8';
                end
                nTrialGroups = trialGroup.count;
                assert(nTrialGroups <= 1);
                if ~nTrialGroups
                    continue
                end
                
                trials = fetch(stimulation.StimTrials*trialGroup,'ORDER BY trial_num');
                conditions = fetch(stimulation.StimConditions*trialGroup,'*');
                
                disp 'Extracting trials...'
                presentations = [];
                for i = 1:length(trials)
                    trial_info = fetch1(stimulation.StimTrials & trials(i), 'trial_params');
                    event = fetch(stimulation.StimTrialEvents & trials(i) & 'event_type="showSubStimulus"','*','ORDER BY event_time');
                    onsets = [event.event_time] + latency;
                    
                    if ~isscalar(onsets)
                        offsets = onsets + median(diff(onsets));
                    else
                        event = fetch(stimulation.StimTrialEvents(trials(i), 'event_type="clearScreen"'),'*');
                        assert(isscalar(event))
                        offsets = event.event_time + latency;
                    end
                    
                    cond = trial_info.conditions;
                    for j=1:length(cond)
                        if onsets(j) > times(1) && offsets(j)+1000 < times(end)
                            s.ori = conditions(cond(j)).condition_info.orientation;
                            s.onset = onsets(j);
                            s.offset = offsets(j);
                            presentations = [presentations; s]; %#ok<AGROW>
                        end
                    end
                end
                durations = double([presentations.offset]-[presentations.onset]);
                assert(std(durations)/mean(durations)<0.03)
                tuple.evoked_bins = floor((mean(durations)+200)/tuple.bin_ms);  % include up to 200 ms offset response
                [~,order] = sort([presentations.onset]);
                presentations = presentations(order);
                nTrials = length(presentations);
                
                % tabulate binned trace segments
                totalTrialDuration = min(diff([presentations.onset]));
                bounds = 0:binDuration:totalTrialDuration;  % bin boundaries
                nBins = length(bounds)-1;
                
                segments = cell(nTrials,1);
                for i=1:nTrials
                    seg = arrayfun(@(j) mean(X(times>presentations(i).onset+bounds(j) & times<presentations(i).onset+bounds(j+1),:)), 1:length(bounds)-1, 'uni', false);
                    segments{i} = cat(1,seg{:});
                end
                [presentations.segment] = deal(segments{:});
                [segments,directions] = dj.struct.tabulate(presentations,'segment','ori');
                nDirs = length(directions);
                segments(cellfun(@isempty, segments)) = {nan(nBins,nCells)};  % replace empties with nans
                tuple.ndirs = nDirs;
                tuple.trace_segments = reshape(single(cat(1,segments{:})), ...
                    nBins, nDirs, [], nCells);
                tuple.directions = directions;
                tuple.ntrials = hist([presentations.ori],directions);
                
                disp inserting..
                self.insert(tuple)
            end
        end
    end
end
