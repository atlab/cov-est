%{
sim.QuadLoss (computed) # quadratic loss
-> sim.CovMatrix
-----
quad_loss : double
%}

classdef QuadLoss < dj.Relvar & dj.AutoPopulate

	properties
		popRel  = sim.CovMatrix & 'nfolds>1'
	end

	methods(Access=protected)

		function makeTuples(self, key)
            
            % X := nBins * nDirs * nTrials * nCells
            [X,R,V] = fetch1(sim.Sample*sim.CovMatrix & key, ...
                'sample_data', 'corr_matrix', 'variances');
            [n,p] = size(X);
            X = double(reshape(X,1,1,n,p));
            
            % split into training and testing
            [~, XTest] = cove.splitTrials(X,key.k,key.nfolds);
            CTest = cove.cov(reshape(XTest,[],p));
            V = diag(sqrt(V(:)));
            key.quad_loss = covest.QuadLoss.loss(V*R*V,CTest);
            self.insert(key)
		end
	end

end