%{
covest.Method (lookup) # covariance matrix estimation methods
method          : tinyint                # covariance estimation method
---
include_spont               : tinyint                       # 1=yes, 0=no
condition                   : tinyint                       # 0=all conditions, 1=first conditoin only, 2=second condition only
regularization              : enum('sample','diag','factor','glasso','lv-glasso') #
loss_fun                    : varchar(255)                  # loss function for model selection
hyperparam_space            : longblob                      # arrays of hyperparameter values
%}

classdef Method < dj.Relvar
    methods
        function fill(self)
            loss =  '@(S,Sigma)(trace(Sigma/S)+covest.lib.logDet(S))/size(S,1)';
            
            self.inserti({
                
            0   1 0  'sample'  loss []
            1   0 0  'sample'  loss []
            2   0 1  'sample'  loss []
            3   0 2  'sample'  loss []
            
            10  1 0  'diag'  loss  {exp(-5:0.05:0),exp(-10:0.05:0)}
            11  0 0  'diag'  loss  {exp(-5:0.05:0),exp(-10:0.05:0)}
            
            20  1 0  'factor'  loss  {0:70,1,0}
            30  1 0  'factor'  loss  {0:70,1,exp(-6:0.05:0)}
            31  0 0  'factor'  loss  {0:70,1,exp(-6:0.05:0)}
            
            40  1 0  'glasso'  loss  {0,exp(-6:0.05:0)}
            41  0 0  'glasso'  loss  {0,exp(-6:0.05:0)}
            
            50  1 0  'glasso'  loss  {exp(-6:.1:-3),exp(-6:0.05:0)}
            
            60  1 0  'lv-glasso' loss {exp(-6:.1:-3),exp(-6:.05:0),exp(-6:.1:1)}
            
            70  1 0  'lv-glasso' loss {0,exp(-6:.05:0),exp(-6:.1:1)}
            71  0 0  'lv-glasso' loss {0,exp(-6:.05:0),exp(-6:.1:1)}
            72  0 1  'lv-glasso' loss {0,exp(-6:.05:0),exp(-6:.1:1)}
            73  0 2  'lv-glasso' loss {0,exp(-6:.05:0),exp(-6:.1:1)}
            
            80  1 0  'glasso'  loss  {0,exp(-6:0.05:0)}
            
            90  1 0  'lv-glasso' loss {0,exp(-6:.05:0),exp(-6:.1:1)}
            91  0 0  'lv-glasso' loss {0,exp(-6:.05:0),exp(-6:.1:1)}
            92  0 1  'lv-glasso' loss {0,exp(-6:.05:0),exp(-6:.1:1)}
            93  0 2  'lv-glasso' loss {0,exp(-6:.05:0),exp(-6:.1:1)}
            })
        end
    end
end
