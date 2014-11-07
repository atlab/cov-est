%{
sim.Rank (computed) # my newest table
-> sim.CovMatrix
-----
rank  : smallint
%}

classdef Rank < dj.Relvar & dj.AutoPopulate

	properties
		popRel  = sim.CovMatrix & 'method=30';
	end

	methods(Access=protected)

		function makeTuples(self, key)
            key.rank = size(fetch1(sim.CovMatrix & key, 'lowrank'),2);
			self.insert(key)
		end
	end

end