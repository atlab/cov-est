%{
covest.CovMatrix (computed) # regularized correlation matrix estimates
-> covest.ActiveCells
-> covest.Method              # method for computing correlation matrices
-> covest.Fold
-----
corr_matrix                  : longblob                      # estimated covariance matrix
variances                    : longblob                      # nBins * nDirs * 1 * nCells
sparse = null                : longblob                      # sparse component of the matrix
lowrank = null               : longblob                      # low-rank component of matrix
hypers = null                : blob                          # values of hyperparameters
visited = null               : longblob                      # tested hyperparamter values
losses  = null               : longblob                      # losses for tested hyperparameter values
sparsity = 0                 : float                         # fraction of zeros off-diagonal
cv_loss = null               : double                        # cross-validation loss
host                         : varchar(255)                  # computer that did the job
cm_ts = CURRENT_TIMESTAMP   : timestamp  
%}

classdef CovMatrix < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel  = covest.ActiveCells * covest.Traces * covest.Method * covest.Fold ...
            & 'ndirs=2 || `condition`=0' & 'corr_estimation="sample"' & 'ncells>50';
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            opt = fetch(covest.Method & key,'*');
            [k,nFolds] = fetch1(covest.Fold & key, 'k','nfolds');
            loss = eval(opt.loss_fun);
            
            % X := nBins * nDirs * nTrials * nCells
            [X,selection] = fetch1(covest.Traces*covest.ActiveCells & key, ...
                'trace_segments', 'selection');
            X = X(:,:,:,selection);
            
            % exclude spontaneous activity if necessary
            if ~opt.include_spont
                evokedBins = fetch1(covest.Traces & key, 'evoked_bins');
                X = X(1:evokedBins,:,:,:);
            end
            
            % select stimulus condition
            if opt.condition
                X = X(:,opt.condition,:,:);
            end
            
            % estimate and subtract mean responses, per condition, per bin
            X = double(X);
            M = nanmean(X,3);
            X = bsxfun(@minus, X, M);
            
            % split into training and testing
            [X, XTest] = covest.splitTrials(X,k,nFolds);
            
            % estimate hyperparameters (if any)
            hypers = {};
            if ~isempty(opt.hyperparam_space)
               [hypers, key.visited, key.losses] = covest.crossEstimateHyper(X, loss, ...
                    opt.var_estimation, opt.corr_estimation, opt.hyperparam_space{:});
            end
            [C,V,extras] = covest.estimate(X, M, opt.var_estimation, opt.corr_estimation, hypers);
            
            key.corr_matrix = C;
            key.variances = V;
            if ~isempty(XTest)
                key.cv_loss = covest.vloss(XTest,C,V,loss);
            end
            if ~isempty(hypers)
                key.hypers = hypers;
            end
            if ~isempty(extras) && ~isempty(fieldnames(extras))
                error 'save extras'
            end
            [~,key.host] = system('hostname');
            self.insert(key)
        end
    end
end