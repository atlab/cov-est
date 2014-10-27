%{
sim.CovMatrix(computed) # estimate of simulated covariance matrix. See also: covest.CovMatrix
-> sim.Sample
-> sim.Method
-> sim.Fold
---
means                       : longblob                      # estimate of the means
variances                   : longblob                      # estimate of variances
corr_matrix                 : longblob                      # estimated covariance matrix
sparse=null                 : longblob                      # sparse component of the matrix
lowrank=null                : longblob                      # low-rank component of matrix
hypers=null                 : blob                          # values of hyperparameters
visited=null                : longblob                      # tested hyperparamter values
losses=null                 : longblob                      # losses for tested hyperparameter values
sparsity=0                  : float                         # fraction of zeros off-diagonal
cv_loss=null                : double                        # cross-validation loss
host                        : varchar(255)                  # computer that did the job
computing_time              : float                         # (s) time required to compute
cm_ts=CURRENT_TIMESTAMP     : timestamp                     #
%}

classdef CovMatrix < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel  =  sim.Sample*sim.Method*sim.Fold
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            t1 = tic;
            opt = fetch(sim.Method & key,'*');
            [k,nFolds] = fetch1(sim.Fold & key, 'k','nfolds');
            
            % X := nBins * nDirs * nTrials * nCells
            X = fetch1(sim.Sample & key, 'sample_data');
            [n,p] = size(X);
            X = double(reshape(X,1,1,n,p));
            
            % split into training and testing
            [X, XTest] = cove.splitTrials(X,k,nFolds);
            
            % estimate hyperparameters (if any)
            if ~any(cellfun(@length,opt.hyperparam_space)>1)
                % hyperparameter values are known
                hypers = cell2mat(opt.hyperparam_space);
            else
                % optimize hyperparameters
                [hypers, ~, key.visited, key.losses] = ...
                    cove.crossEstimateHyper(X, 1, ...
                    opt.regularization, opt.hyperparam_space);
            end
            [R,M,V,extras] = cove.estimate(X, 1, opt.regularization, hypers);
            key.corr_matrix = R;
            key.means = M;
            key.variances = V;
            if ~isempty(hypers)
                key.hypers = hypers;
            end
            if ~isempty(XTest)
                key.cv_loss = cove.vloss(XTest, R, M, V, 0);
            end
            if ~isempty(extras) && ~isempty(fieldnames(extras))
                if isfield(extras,'loading_matrix')
                    key.lowrank = extras.loading_matrix;
                end
                if isfield(extras,'H') && numel(extras.H)>0
                    key.lowrank = extras.H;
                end
                if isfield(extras,'S')
                    key.sparsity = cove.sparsity(extras.S);
                    key.sparse = extras.S;
                end
            end
            [~,key.host] = system('hostname');
            key.computing_time = toc(t1);
            self.insert(key)
        end
    end
    
end