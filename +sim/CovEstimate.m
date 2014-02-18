%{
sim.CovEstimate (computed) # covariance matrices estimated on 10 folds
-> sim.Sample
-> sim.CovEstimator
-----
cov_matrix  : longblob   #  covariance matrix estimate
excess_loss : float      # loss compared to real truth, only in simulation
hypers      : blob       # values of hyperparameters
%}

classdef CovEstimate < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = sim.Sample*sim.CovEstimator
    end
    
    properties(Constant)
        loss = @(S,Sigma) (trace(Sigma/S) + logDet(S))/size(S,1);
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            Sigma = fetch1(sim.TrueCov & key, 'covmatrix');
            X = sim.Sample.makeSample(key.seed, key.sample_size, Sigma);
            [estimator,searchSpace] = fetch1(sim.CovEstimator & key,...
                'estimator_call','hyperparam_space');
            if isempty(searchSpace)
                C = covest.estimate(X, estimator);
                hypers = [];
            else
                [C, hypers] = ...
                    covest.crossEstimateHyper(X, self.loss, estimator, searchSpace{:});
            end
            key.cov_matrix = single(C);
            key.excess_loss = self.loss(C,Sigma)-self.loss(Sigma,Sigma);
            key.hypers = hypers;
            self.insert(key)
        end
    end
end
