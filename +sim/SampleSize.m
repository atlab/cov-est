%{
sim.SampleSize (lookup) # sample size
sample_size  :  smallint  #  sample size
-----

%}

classdef SampleSize < dj.Relvar
    
    properties(Constant)
        table = dj.Table('sim.SampleSize')
    end
    
    methods
        function fill(self)
            for n = [150 250 500 1000 2000 4000]
                self.inserti(struct('sample_size',n))
            end
        end
    end
end