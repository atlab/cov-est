%{
sim.Sample (computed) # random sample
-> sim.Truth
-> sim.SampleSize
-> sim.Model
sample_seed     : smallint               # 
---
sample_data                 : longblob                      # raw data
%}

classdef Sample < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = sim.Truth*sim.SampleSize*sim.Model
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            model = fetch1(sim.Model & key, 'model_name');            
            % make sure the model seed differs from the truth seed
            sampleSeed = 32000+fetch1(sim.Truth & key, 'truth_seed');
            assert(sampleSeed ~= fetch1(sim.Truth & key, 'truth_seed'))
            rng(sampleSeed)
            n = fetch1(sim.SampleSize & key, 'sample_size');
            C = fetch1(sim.Truth & key, 'true_cov');
            key.sample_seed = sampleSeed;
            p = size(C,1);
            switch model
                case 'gauss'
                    key.sample_data = single(mvnrnd(zeros(p,1),C,n));
                case 'ising'
                    key.sample_data = int8(ising(C,zeros(p,1),n));
            end
            self.insert(key)
        end
    end
    
end



function X = ising(C,h,n)
% sample from Ising model using Metropolis-Hastings

p = size(C,1);
S = -inv(C);
S = S - diag(diag(S));
x = ones(1,p);
skip = 1000;   % it's kinda arbitrary
X = nan(n,p);
for i=1:n
    for t=1:skip
        xnew = randi(2,1,p)*2-3;   % this works well for small problems
        probRatio = exp(xnew*S*xnew' - x*S*x');
        if rand < probRatio
            x = xnew;
        end
    end
    X(i,:) = x;
end
end