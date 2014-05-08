%{
sim.SampleSize (lookup) # sample size (recording duration)
sample_size  : smallint   # sample size (recording duration)
-----
%}

classdef SampleSize < dj.Relvar
    methods
        function fill(self)
            self.insert({
                250
                500
                1000
                2000
                4000
                })
        end
    end
end