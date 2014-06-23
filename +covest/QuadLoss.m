%{
covest.QuadLoss (computed) # my newest table
-> covest.CovMatrix
-----
quad_loss=null  :  double   #  quadratic loss value
%}

classdef QuadLoss < dj.Relvar & dj.AutoPopulate
    
    properties(Constant)
        loss = @(S,Sigma) trace(Sigma/S-eye(size(S,1)))/size(S,1)^2;
    end
    
    properties
        popRel  = covest.CovMatrix & 'nfolds>1'
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            opt = fetch(covest.Method & key,'*');
            [k,nFolds] = fetch1(covest.Fold & key, 'k','nfolds');
            
            % X := nBins * nDirs * nTrials * nCells
            [X,selection] = fetch1(covest.Traces*covest.ActiveCells & key, ...
                'trace_segments', 'selection');
            X = double(X(:,:,:,selection));
            [nBins, nDirs, nTrials, nCells] = size(X);
            
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
            [M,V,R,delta] = fetch1(covest.CovMatrix & key, ...
                'means','variances','corr_matrix','delta');
            [~, XTest] = cove.splitTrials(X,k,nFolds);
            assert(~isempty(XTest))
            X = bsxfun(@minus, X, M);            
            V0 = reshape(nanmean(reshape(V,[],nCells)),[1 1 1 nCells]);
            V = bsxfun(@plus,(1-delta)*V,delta*V0);
            X = bsxfun(@rdivide, X, sqrt(V));
            RTest = cov(reshape(X,[],nCells));
            
            
            
            key.quad_loss = self.loss(R,RTest);
            self.insert(key)
        end
    end
    
end