%{
sim.CovEstimateFold (computed) # covariance matrices estimated on 10 folds
-> sim.Sample
-> sim.CovEstimator
-> sim.Fold
-----
cov_matrix  : longblob   # covariance matrix estimate
cross_loss  : float      # loss estimated from validation subset
hypers      : blob       # values of hyperparameters
%}

classdef CovEstimateFold < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = sim.Sample*sim.CovEstimator*sim.Fold
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            Sigma = fetch1(sim.TrueCov & key, 'covmatrix');
            X = sim.Sample.makeSample(key.seed, key.sample_size, Sigma);
            [n,p] = size(X);
            K = 10;   % cross-validation fold
            rng(0)
            cvInd = covest.crossvalIndex(n,K,round(n/(3*K)));
            XTrain = X(cvInd~=key.fold,:);
            CTest  = covest.cov(X(cvInd==key.fold,:));  % testing set
            
            loss = @(S,Sigma) (trace(Sigma/S) + logDet(S))/size(S,1);
            [estimator,searchSpace] = fetch1(sim.CovEstimator & key,...
                'estimator_call','hyperparam_space');
            if isempty(searchSpace)
                C = covest.estimate(XTrain,  estimator);
                hypers = [];
            else
                [C, hypers] = ...
                    covest.crossEstimateHyper(XTrain, sim.CovEstimate.loss, estimator, searchSpace{:});
            end
            key.cov_matrix = C;
            key.cross_loss = loss(C,CTest);
            key.hypers = hypers;
            
            self.insert(key)
        end
    end
end
