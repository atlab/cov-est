%{
pop.AodCovEstimateFold (computed) # covariance matrix estimate
-> pop.AodBinnedTraces
-> pop.CovEstimator
-> pop.Fold
---
covmatrix                   : longblob                      # estimated covariance matrix
cross_loss=null             : float                         # value of the loss function
hypers                      : blob                          # optimal values of hyperparameters
%}

classdef AodCovEstimateFold < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = pop.AodBinnedTraces*pop.CovEstimator*pop.Fold
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            X = fetch1(pop.AodBinnedTraces & key, 'binned_traces');
            X = double(X);
            [n,p] = size(X);
            
            [estimator,searchSpace,loss] = fetch1(pop.CovEstimator & key,...
                'estimator_call','hyperparam_space','loss_fun');
            loss = eval(loss);
            [k,K] = fetch1(pop.Fold & key, 'k','nfolds');
            
            rng(0)
            cvInd = covest.crossvalIndex(n,K,round(n/(3*K)));
            XTrain = X(cvInd~=k,:);  % training set
            CTest  = covest.cov(X(cvInd==k,:));  % testing covariance
            
            disp 'cross-validation'
            if isempty(searchSpace)
                C = covest.estimate(XTrain,  estimator);
                hypers = [];
            else
                [C, hypers] = ...
                    covest.crossEstimateHyper(XTrain, sim.CovEstimate.loss, estimator, searchSpace{:});
            end
            key.covmatrix = C;
            key.cross_loss = loss(C,CTest);
            key.hypers = hypers;
            
            self.insert(key)
        end
    end
end
