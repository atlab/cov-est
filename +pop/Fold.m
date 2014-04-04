%{
pop.Fold (lookup) # my newest table
k      :   tinyint   # kth fold
nfolds :   tinyint   # total folds
-----

%}

classdef Fold < dj.Relvar

	properties(Constant)
		table = dj.Table('pop.Fold')
    end
    
    methods
        function fill(self)
            n = 10;
            for k=1:n
                self.inserti(struct('k',k,'nfolds',n))
            end
        end
    end
end
