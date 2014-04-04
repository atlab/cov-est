%{
pop.AodCovEstimate (computed) # covariance matrix estimate
-> pop.AodBinnedTraces
-> pop.CovEstimator
---
cov_matrix                   : longblob                      # estimated covariance matrix
sparse = null                : longblob                      # sparse component of the matrix
lowrank = null               : longblob                      # low-ranke component of matrix
hypers = null                : blob                          # values of hyperparameters
sparsity = 0                 : float                         # fraction of zeros off-diagonal
%}

classdef AodCovEstimate < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = pop.AodBinnedTraces*pop.CovEstimator & 'bin_opt=0'
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            X = fetch1(pop.AodBinnedTraces & key, 'binned_traces');
            X = double(X);
            rng(0);  % make it fair
            
            [estimator,searchSpace,loss] = fetch1(pop.CovEstimator & key,...
                'estimator_call','hyperparam_space','loss_fun');
            loss = eval(loss);
            
            disp 'cross-validation'
            if isempty(searchSpace)
                C = covest.estimate(X,  estimator);
            else
                [C, hypers, extras] = ...
                    covest.crossEstimateHyper(X, loss, estimator, searchSpace{:});
            end
            
            key.cov_matrix = single(C);
            if ~isempty(whos('hypers'))
                key.hypers = hypers;
            end
            if ~isempty(whos('extras'))
                if isfield(extras,'S')
                    key.sparse = single(extras.S);
                    key.sparsity = covest.sparsity(extras.S);
                end
                if isfield(extras,'loading_matrix')
                    key.lowrank = single(extras.loading_matrix);
                elseif isfield(extras,'H')
                    key.lowrank = single(extras.H);
                end
            end
            self.insert(key)
        end
    end
end