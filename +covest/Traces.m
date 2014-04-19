% The contents of this table are copied from computed table pop.TrialTraceSet

%{
covest.Traces (manual) # traces associated with trials, aligned on stimulus onset
subject_id      : int unsigned           # unique identifier for subject
setup           : tinyint unsigned       # setup number
session_start_time: bigint               # start session timestamp
aod_scan_start_time: bigint              # start session timestamp
preprocess_method_num: int unsigned      # Preprocessing method
high_repeats                : tinyint                       # 1=this dataset has few stimulus conditions with many trials
---
latency_ms                  : smallint                      # (ms) presumed screen-to-V1 latency
bin_ms                      : float                         #
evoked_bins                 : tinyint                       # number of evoked bins -- the remainder are spontaneous
nneurons                    : smallint                      # the number of neurons included in the analysis
cell_xyz                    : blob                          #
cellnums                    : blob                          # cell selection
ndirs                       : smallint                      # number of directions
directions                  : blob                          # directions
ntrials                     : blob                          # trial in each direction
trace_segments              : longblob                      # trace segements: nBins x nDirs x nTrials x nCells
%}

classdef Traces < dj.Relvar
end