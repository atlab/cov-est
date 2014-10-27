%{
covest.CovMatrix (computed) # regularized correlation matrix estimates
-> covest.ActiveCells
-> covest.Method
-> covest.Fold
---
means                       : longblob                      # estimate of the means
variances                   : longblob                      # estimate of variances
corr_matrix                 : longblob                      # estimated covariance matrix
test_matrix=null            : longblob                      # sample cov matrix from the testing set
sparse=null                 : longblob                      # sparse component of the matrix
lowrank=null                : longblob                      # low-rank component of matrix
hypers=null                 : blob                          # values of hyperparameters
delta                       : float                         # variance regularization
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
        popRel  = covest.ActiveCells * covest.Traces * covest.Method * covest.Fold ...
            & 'preprocess_method_num=5' & 'ndirs=2 || `condition`=0' & 'ncells>100'
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            t1 = tic;
            opt = fetch(covest.Method & key,'*');
            [k,nFolds] = fetch1(covest.Fold & key, 'k','nfolds');
            
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
            [X, XTest] = cove.splitTrials(X,k,nFolds);
            
            if any(cellfun(@length,opt.hyperparam_space)>1)
                % optimize hyper parameters
                [hypers, bestDelta, key.visited, key.losses] = ...
                    cove.crossEstimateHyper(X, evokedBins, ...
                    opt.regularization, opt.hyperparam_space);
            else
                % no hyperparameters. Optimize delta (variance regularizer)
                hypers = cell2mat(opt.hyperparam_space);
                K = 10;  % inner-loop cross-validation
                [XTest_,R_,M_,V_] = arrayfun(@(k) estimate_(hypers,k,K), 1:K, 'uni', false);
                bestDelta = mean(cellfun(@(xt,R,M,V) ...
                    cove.findBestDelta(xt,R,M,V), ...
                    XTest_, R_, M_, V_));
            end
            [R,M,V,extras] = cove.estimate(X, evokedBins, opt.regularization, hypers);
            key.corr_matrix = R;
            key.means = M;
            key.variances = V;
            key.delta = bestDelta;
            if ~isempty(hypers)
                key.hypers = hypers;
            end
            if ~isempty(XTest)
                key.test_matrix = XTest; 
                key.cv_loss = cove.vloss(XTest, R, M, V, bestDelta);
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
            
            
            function [XTest,R,M,V] = estimate_(hypers,k,K)
                % compute cross-validation loss
                [XTrain,XTest] = cove.splitTrials(X,k,K);
                [R, M, V] = cove.estimate(XTrain, evokedBins, opt.regularization, hypers);
            end
        end
    end
end
