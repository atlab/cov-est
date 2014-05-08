%{
sim.Method (lookup) # covariance matrix estimation methods
method          : tinyint                # covariance estimation method
---
regularization              : enum('sample','diag','factor','glasso','lv-glasso') #
loss_fun                    : varchar(255)                  # loss function for model selection
hyperparam_space            : longblob                      # arrays of hyperparameter values
%}

classdef Method < dj.Relvar
    methods
        function fill(self)
            loss =  '@(S,Sigma)(trace(Sigma/S)+cove.logDet(S))/size(S,1)';
            
            self.inserti({
                
            0   'sample'    loss  []
            10  'diag'      loss  {exp(-5:0.05:0),exp(-8:0.05:0)}
            30  'factor'    loss  {0:20,1,exp(-6:0.02:0)}
            80  'glasso'    loss  {0,exp(-6:0.05:0)}
            90  'lv-glasso' loss  {0,exp(-6:.05:0),exp(-6:.1:1)}
            
            })
        end
    end
end
