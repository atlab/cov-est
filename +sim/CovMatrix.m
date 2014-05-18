%{
sim.CovMatrix (computed) # estimate of simulated covariance matrix. See also: covest.CovMatrix
-> sim.Sample
-> sim.Method
-> covest.Fold
---
cov_matrix                  : longblob                      # estimated covariance matrix
test_matrix=null            : longblob                      # sample cov matrix from the testing set
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
        popRel  =  sim.Sample*sim.Method*covest.Fold
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            t1 = tic;
            opt = fetch(sim.Method & key,'*');
            [k,nFolds] = fetch1(covest.Fold & key, 'k','nfolds');
            loss = eval(opt.loss_fun);
            
            % X := nBins * nDirs * nTrials * nCells
            X = fetch1(sim.Sample & key, 'sample_data');
            [n,p] = size(X);
            X = double(reshape(X,1,1,n,p));
            
            % split into training and testing
            [X, XTest] = cove.splitTrials(X,k,nFolds);
            
            % estimate hyperparameters (if any)
            hypers = {};
            if ~isempty(opt.hyperparam_space)
                [hypers, key.visited, key.losses] = cove.crossEstimateHyper(X, 1, loss, ...
                    opt.regularization, opt.hyperparam_space);
            end
            [C,~,extras] = cove.estimate(X, 0, 1, opt.regularization, hypers);
            
            key.cov_matrix = C;
            if ~isempty(XTest)
                CTest = cove.estimate(XTest, 0, 1, 'sample', {});
                key.test_matrix = CTest;
                key.cv_loss = loss(C,CTest);
            end
            if ~isempty(hypers)
                key.hypers = hypers;
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