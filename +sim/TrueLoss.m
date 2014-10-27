%{
sim.TrueLoss (computed) # true loss based on ground truth
-> sim.CovMatrix
-----
true_loss  : longblob
%}

classdef TrueLoss < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = sim.CovMatrix
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            loss = @(C,trueC) (trace(trueC/C)+cove.logDet(C))/size(C,1);
            [trueC,R,V] = fetch1(pro(sim.Truth,'true_cov')*sim.CovMatrix & key, ...
                'true_cov', 'corr_matrix','variances');
            V = diag(sqrt(V(:)));
            key.true_loss = loss(V*R*V,trueC)-loss(trueC,trueC);
            self.insert(key)
        end
    end
    
end 