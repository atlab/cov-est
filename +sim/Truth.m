%{
sim.Truth (manual) # ground truth covariance matrices
truth           : smallint
---
population_size             : smallint                      # number of neurons
truth_seed                  : smallint                      # random number generate seed
truth_type                  : enum('diag','factor','sparse','sparse+latent') #
sparsity                    : double                        # sparsity of the sparse component 
nfactors                    : smallint                      # rank of the low-rank component 
true_cov                    : longblob                      # true covariance matrix
%}

classdef Truth < dj.Relvar
    methods
        function fill(self)
            p = 50;  % population size
            id = 0;
            for truthType = {'diag','factor','sparse','sparse+latent'}
                for seed = 1:30
                    % generate the seed matrix
                    rng(seed)
                    V =  lognrnd(1,.25,50,1);  % variances
                    X = mvnrnd(zeros(size(V))',diag(V),p*3);
                    C = cove.cov(X);
                    nFactors = p;
                    sparsity = 0;
                    
                    switch truthType{1}
                        case 'diag'
                            C = diag(diag(C));
                        
                        case 'factor'
                            nFactors = 4;
                            [L,V] = cove.factor(C,nFactors);
                            C = L*L' + diag(V);
                            
                        case 'sparse'
                            nFactors = 0;
                            res = cove.lvglasso(C,0.4,0.2,struct('refit',true,'max_latent',0));
                            sparsity = cove.sparsity(res.S);
                            C = inv(res.S);

                        case 'sparse+latent'
                            nFactors = 3;
                            res = cove.lvglasso(C,0.4,0.1,struct('refit',true,'max_latent',nFactors));
                            sparsity = cove.sparsity(res.S);
                            nFactors = size(res.eigL,1);
                            C = inv(res.S-res.L);
                            
                        otherwise 
                            error 'unknown truth type'
                    end
                    
                    id = id + 1;
                    self.inserti({
                        id  p seed truthType{1} sparsity nFactors C
                        })
                end
            end
        end
    end
end
