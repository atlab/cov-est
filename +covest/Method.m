%{
covest.Method (lookup) # covariance matrix estimation methods
method          : tinyint                # covariance estimation method
---
include_spont             : tinyint                       # 1=yes, 0=no
condition                 : tinyint                       # 0=all conditions, 1=first conditoin only, 2=second condition only
var_regularization        : enum('','L2')                 # regularization of variances
cov_regularization        : enum('sample','diag','factor','glasso','lv-glasso') #
hyperparam_space          : longblob                      # arrays of hyperparameter values
%}

classdef Method < dj.Relvar
    methods
        function fill(self)
            
            self.inserti({
                
            % common-variance sample covariances
            0   1 0  '' 'sample'  {}
            1   0 0  '' 'sample'  {}
            2   0 1  '' 'sample'  {}
            3   0 2  '' 'sample'  {}
            
            % regularized variance and covariance
            5   1 0  'L2' 'sample'  {exp(-1.5:0.03:0)}
            
            10  1 0   ''  'diag'                     {exp(-6:0.05:-1) exp(-4:0.05:0)}
            15  1 0  'L2' 'diag'    {exp(-1.5:0.03:0) exp(-6:0.05:-1) exp(-4:0.05:0)}
            
            30  1 0   ''  'factor'                   {exp(-6:0.05:-1) 0:70}
            35  1 0  'L2' 'factor'  {exp(-1.5:0.03:0) exp(-6:0.05:-1) 0:70}
            
            80  1 0   ''  'glasso'                   {exp(-5:0.05:-1)}
            85  1 0  'L2' 'glasso'  {exp(-1.5:0.04:0) exp(-5:0.05:-1)}
            
            90  1 0   ''  'lv-glasso'                  {exp(-5:.05:-1) exp(-4:.05:-1)}
            95  1 0  'L2' 'lv-glasso' {exp(-1.5:0.03:0) exp(-5:.05:-1) exp(-4:.05:-1)}
            
            })
        end
    end
end
