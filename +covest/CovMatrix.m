%{
covest.CovMatrix (computed) # regularized correlation matrix estimates
-> covest.ActiveCells
-> covest.Method              # method for computing correlation matrices
-> covest.Fold
-----
cov_matrix                   : longblob                      # estimated covariance matrix
sparse = null                : longblob                      # sparse component of the matrix
lowrank = null               : longblob                      # low-rank component of matrix
hypers = null                : blob                          # values of hyperparameters
visited = null               : longblob                      # tested hyperparamter values
losses  = null               : longblob                      # losses for tested hyperparameter values
sparsity = 0                 : float                         # fraction of zeros off-diagonal
cv_loss = null               : double                        # cross-validation loss
host                         : varchar(255)                  # computer that did the job
computing_time               : float                         # (s) time required to compute
cm_ts = CURRENT_TIMESTAMP   : timestamp
%}

classdef CovMatrix < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel  = covest.ActiveCells * covest.Traces * covest.Method * covest.Fold ...
            & 'preprocess_method_num=5' & 'ndirs=2 || `condition`=0' & 'ncells>100'
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            t1 = tic;
            opt = fetch(covest.Method & key,'*');
            [k,nFolds] = fetch1(covest.Fold & key, 'k','nfolds');
            loss = eval(opt.loss_fun);
            
            % X := nBins * nDirs * nTrials * nCells
            [X,selection] = fetch1(covest.Traces*covest.ActiveCells & key, ...
                'trace_segments', 'selection');
            X = double(X(:,:,:,selection));
            
            % exclude spontaneous activity if necessary
            evokedBins = fetch1(covest.Traces & key, 'evoked_bins');
            if ~opt.include_spont
                X = X(1:min(evokedBins,end),:,:,:);
            end
            
            % select stimulus condition
            if opt.condition
                X = X(:,opt.condition,:,:);
            end
            
            % split into training and testing
            [X, XTest] = covest.splitTrials(X,k,nFolds);
            
            % estimate hyperparameters (if any)
            hypers = {};
            if ~isempty(opt.hyperparam_space)
                [hypers, key.visited, key.losses] = covest.crossEstimateHyper(X, evokedBins, loss, ...
                    opt.regularization, opt.hyperparam_space);
            end
            [C,~,extras] = covest.estimate(X, [], evokedBins, opt.regularization, hypers);
            
            key.cov_matrix = C;
            if ~isempty(XTest)
                CTest = covest.estimate(XTest,[],evokedBins, 'sample', {});
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
                    key.sparsity = covest.lib.sparsity(extras.S);
                    key.sparse = extras.S;
                end
            end
            [~,key.host] = system('hostname');
            key.computing_time = toc(t1);
            self.insert(key)
        end
    end
end