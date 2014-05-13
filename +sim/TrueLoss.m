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
            loss = @(C,trueC) (trace(C/trueC)+cove.logDet(trueC))/size(trueC,1);
            [trueC,C] = fetch1(pro(sim.Truth,'true_cov')*sim.CovMatrix & key, 'true_cov', 'cov_matrix');
            key.true_loss = loss(trueC,C)-loss(trueC,trueC);
            self.insert(key)
        end
    end
    
end