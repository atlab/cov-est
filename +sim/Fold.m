%{
sim.Fold (lookup) # helps parallelize cross-validation
fold  : tinyint # cross-validation fold number, 0 - entire dataset, 1-10 folds
-----
%}

classdef Fold < dj.Relvar
    
    properties(Constant)
        table = dj.Table('sim.Fold')
    end
    
    methods
        function fill(self)
            for fold=0:10
                self.insert(struct('fold',fold))
            end
        end
    end
end