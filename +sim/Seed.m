%{
sim.Seed (lookup) # my newest table
seed : tinyint # rng seed
%}

classdef Seed < dj.Relvar

	properties(Constant)
		table = dj.Table('sim.Seed')
    end
    
    
    methods
        function fill(self)
            for seed = 1:30
                self.inserti(struct('seed',seed))
            end
        end
    end
end
