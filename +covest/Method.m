%{
covest.Method (lookup) # covariance matrix estimation methods
method          : tinyint                # covariance estimation method
---
include_spont             : tinyint                       # 1=yes, 0=no
condition                 : tinyint                       # 0=all conditions, 1=first conditoin only, 2=second condition only
mean_regularization       : enum('','L1','L2')               # regularization of the mean response
cov_regularization        : enum('sample','diag','factor','glasso','lv-glasso') # 
hyperparam_space          : longblob                      # arrays of hyperparameter values
%}

classdef Method < dj.Relvar
    methods
        function fill(self)

            self.inserti({
                
            % common-variance sample covariances
            0   1 0  '' 'sample'  {1}
            1   0 0  '' 'sample'  {1}
            2   0 1  '' 'sample'  {1}
            3   0 2  '' 'sample'  {1}

            % regularized means and covariance
            4   1 0  'L1' 'sample'  {exp(-5:0.05:0) exp(-5:0.05:0)}
            5   1 0  'L2' 'sample'  {exp(-5:0.05:0) exp(-5:0.05:0)}
            
            10  1 0  'L1' 'diag'  {exp(-5:0.05:0) exp(-5:0.05:0) exp(-5:0.05:0) exp(-10:0.05:0)}
            15  1 0  'L2' 'diag'  {exp(-5:0.05:0) exp(-5:0.05:0) exp(-5:0.05:0) exp(-10:0.05:0)}
            
%             30  1 0  'L2' 'factor'  {0:70 1 exp(-6:0.05:0)}
%             
%             80  1 0  'L2' 'glasso'  {0 exp(-6:0.05:0)}
% 
%             90  1 0  'L2' 'lv-glasso' {0 exp(-6:.05:0) exp(-6:.1:1)}
            
            })
        end
    end
end
