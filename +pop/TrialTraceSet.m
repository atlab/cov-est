%{
pop.TrialTraceSet (computed) # traces associated with trials
-> pop.StableScans
-> aod.UniqueCells
-> aod.TracePreprocessSet
-> pop.BinOpt
-----
nneurons                    : smallint                      # the number of neurons included in the analysis
cellnums                    : longblob                      # cell selection
ndirs                       : smallint                      # number of directions
dirs                        : longblob                      # direction of grating movement
%}

classdef TrialTraceSet < dj.Relvar & dj.AutoPopulate

	properties
        popRel = (aod.TracePreprocessSet*aod.UniqueCells*pop.StableScans*pop.BinOpt) ...
            & acq.AodStimulationLink & 'preprocess_method_num=5' & 'bin_opt=0'
	end

	methods(Access=protected)

		function makeTuples(self, key)
			self.insert(key)
		end
	end

end