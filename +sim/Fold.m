%{
sim.Fold (lookup) # my newest table
k      :   tinyint   # kth fold
nfolds :   tinyint   # total folds
-----
%}

classdef Fold < dj.Relvar
     
    methods
        function fill(self)
            self.inserti({1 1});
            n = 10;
            for k=1:n
                self.inserti({k n})
            end
        end
    end
end