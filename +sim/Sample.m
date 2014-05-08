%{
sim.Sample (computed) # random sample
-> sim.Truth
-> sim.SampleSize
sample_seed  : tinyint
-----
sample_data  :  longblob    # raw data
%}

classdef Sample < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = sim.Truth*sim.SampleSize
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            sampleSeed = 1;
            rng(sampleSeed)
            n = fetch1(sim.SampleSize & key, 'sample_size');
            C = fetch1(sim.Truth & key, 'true_cov');
            key.sample_seed = sampleSeed;
            key.sample_data = single(mvnrnd(zeros(size(C,1),1),C,n));
            self.insert(key)
        end
    end
    
end