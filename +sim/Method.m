%{
sim.Method (lookup) # covariance matrix estimation methods
method          : tinyint                # covariance estimation method
---
regularization              : enum('sample','diag','factor','glasso','lv-glasso') #
hyperparam_space            : longblob                      # arrays of hyperparameter values
%}

classdef Method < dj.Relvar
    methods
        function fill(self)
            
            self.inserti({
                
            0   'sample'    {}
            10  'diag'      {exp(-5:0.05:0),exp(-5:0.05:0)}
            30  'factor'    {exp(-6:0.05:-1) 0:30}
            80  'glasso'    {exp(-5:0.05:-1)}
            90  'lv-glasso' {exp(-7:.05:-1.5) exp(-5:.05:-1.5)}
            
            })
        end
    end
end