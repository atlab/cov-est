%{
pop.TrialTraceSet (computed) # traces associated with trials
-> pop.StableScans
-> aod.UniqueCells
-> aod.TracePreprocessSet
-> pop.BinOpt
-----
nneurons                    : smallint                      # the number of neurons included in the analysis
cellnums                    : longblob                      # list of participating cell numbers
ndirs=null                  : smallint                      # number of directions
dirs=null                   : longblob                      # direction of grating movement
psths=null                  : longblob                      # peristimulus histograms for each cell and direction
binned_traces               : longblob                      # array of binned traces
binsize_ms                  : float                         # true bin size
epoch_seconds               : float                         # the duration of the analyzed episode

%}

classdef TrialTraceSet < dj.Relvar & dj.AutoPopulate

	properties
		popRel  % !!! update the populate relation
	end

	methods(Access=protected)

		function makeTuples(self, key)
		%!!! compute missing fields for key here
			self.insert(key)
		end
	end

end