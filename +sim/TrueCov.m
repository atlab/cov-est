%{
sim.TrueCov (lookup) # true covariance matrices used in simulations
true_cov_id  : smallint   # true covariance matrix 
-----
p            : smallint  # population size
covclass     : enum('unstructured','independent','multifactor','sparse','sparse+lowrank') # the structure of the true correlation matrix
covmatrix    : longblob  # the true covariance matrix
nlatent = 0  : smallint  # number of latent models in factor and lv-glasso
sparsity = 0 : float     # fraction of
stddev       : float     # standard deviation of correlations
%}

classdef TrueCov < dj.Relvar
    
    methods
        function fill(self)
            %self.fill_(50,6,0);
            self.fill_(50,20,2);
            self.fill_(50,30,3);            
        end
        
        function id = fill_(self,p,id,seed)
            for name = {'unstructured','independent','multifactor','sparse','sparse+lowrank'}
                id = id+1;                
                disp(name{1})
                sparsity = 0;
                nlatent = 0;
                rng(seed)
                n = 150;  % number of samples in the prematrix -- regulates scale of correlations
                    Sigma = covest.cov(mvnrnd(zeros(p,1),eye(p),n));
                switch name{1}
                    case 'unstructured'
                        % do nothing
                         
                    case 'independent'
                        Sigma = diag(diag(Sigma));
                        
                    case 'multifactor'
                        nlatent = 2;
                        [L,Ph] = covest.factor(Sigma,nlatent);
                        Sigma = L*L' + diag(Ph);
                        
                    case 'sparse'
                        ret = covest.lvglasso(Sigma, 0.12, 0.0, struct('refit',true,'mu',p,'max_latent',0));
                        Sigma = inv(ret.S-ret.L);
                        sparsity = covest.sparsity(ret.S);
                        
                    case 'sparse+lowrank'
                        nlatent = 1;
                        ret = covest.lvglasso(Sigma, 0.12, 0.0, struct('refit',true,'mu',p,'max_latent',nlatent));
                        Sigma = inv(ret.S-ret.L);
                        sparsity = covest.sparsity(ret.S);
                    otherwise
                        error 'something''s wrong'
                end
                imagesc(Sigma,[-1 1]*.25), colormap(covest.doppler); axis image; colorbar
                C = corrcov(Sigma);
                [x,y] = meshgrid(1:p,1:p);
                self.inserti(struct(...
                    'true_cov_id',id,...
                    'covclass', name{1},...
                    'p', p, ...
                    'covmatrix', Sigma,...
                    'nlatent', nlatent,...
                    'sparsity', sparsity,...
                    'stddev', std(C(x<y))));
            end
        end
    end
end
