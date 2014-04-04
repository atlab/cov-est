%{
pop.CVLoss (computed) # cross-validation loss$
-> pop.AodCovEstimateFold2
---
quadratic_loss              : double                        # quadratic covariance matrix loss
normal_loss=null            : double                        # multivariate normal loss
%}

classdef CVLoss < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = pop.AodCovEstimateFold
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            X = fetch1(pop.AodBinnedTraces & key, 'binned_traces');
            X = double(X);
            [n,p] = size(X);
            [k,K] = fetch1(pop.Fold & key, 'k','nfolds');
            
            rng(0)
            cvInd = covest.crossvalIndex(n,K,round(n/(3*K)));
            CTest  = covest.cov(X(cvInd==k,:));  % testing covariance
            
            C = fetch1(pop.AodCovEstimateFold & key, 'covmatrix');
            key.normal_loss = (trace(CTest/C)+logDet(C))/p;
            key.quadratic_loss = trace(CTest/C-eye(p))^2/p^2;
            
            self.insert(key)
        end
    end
end
