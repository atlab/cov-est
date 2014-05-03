%{
covest.QuadLoss (computed) # my newest table
-> covest.CovMatrix
-----
quandratic_loss  :  double   #  quadratic loss value
%}

classdef QuadLoss < dj.Relvar & dj.AutoPopulate

	properties
		popRel  = covest.CovMatrix
	end

	methods(Access=protected)

		function makeTuples(self, key)
		%!!! compute missing fields for key here
			self.insert(key)
		end
	end

end