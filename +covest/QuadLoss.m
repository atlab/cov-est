%{
covest.QuadLoss (computed) # my newest table
-> covest.CovMatrix
-----
quad_loss  :  double   #  quadratic loss value
%}

classdef QuadLoss < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel  = covest.CovMatrix & 'nfolds>1'
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            opt = fetch(covest.Method & key,'*');
            [k,nFolds] = fetch1(covest.Fold & key, 'k','nfolds');
            loss = @(S,Sigma) trace(Sigma/S-eye(size(S,1)))/size(S,1)^2;
            
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
            C = fetch1(covest.CovMatrix & key, 'cov_matrix');
            [~, XTest] = covest.splitTrials(X,k,nFolds);
            CTest = covest.estimate(XTest,[],evokedBins, 'sample', {});
            key.quad_loss = loss(C,CTest);
            self.insert(key)
        end
    end
    
end