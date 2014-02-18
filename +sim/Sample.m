%{
sim.Sample (computed) # my newest table
-> sim.TrueCov
-> sim.SampleSize
-> sim.Seed
---
sample_covmatrix            : longblob                      # sample covariance matrix
%}

classdef Sample < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = sim.TrueCov * sim.SampleSize * sim.Seed
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            Sigma = fetch1(sim.TrueCov & key, 'covmatrix');
            X = sim.Sample.makeSample(key.seed, key.sample_size, Sigma);
            key.sample_covmatrix = covest.cov(X);
            self.insert(key)
        end
    end
    
    methods(Static)
        function X = makeSample(seed,sampleSize,Sigma)
            rng(seed)
            X = mvnrnd(zeros(size(Sigma,1),1),Sigma,sampleSize);
        end
    end
end
