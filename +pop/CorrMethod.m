%{
pop.CorrMethod (lookup) # correlation matrix estimation methods
corr_method     : tinyint                # covarianc
---
include_spont               : tinyint                       # 1=yes, 0=no
condition                   : tinyint                       # 0=all conditions, 1=first conditoin only, 2=second condition only
var_estimation              : enum('uniform','linear to mean','per condition','per bin') # ways to estimate the variance
corr_estimation             : enum('sample','diag','factor','lv-glasso') #
loss_fun                    : varchar(255)                  # loss function for model selection
hyperparam_space            : longblob                      # arrays of hyperparameter values
%}

classdef CorrMethod < dj.Relvar
    methods
        function fill(self)
            loss =  '@(S,Sigma)(trace(Sigma/S)+covest.logDet(S))/size(S,1)';
            self.inserti({
                
            0   0 0  'uniform'        'sample'  loss []
            1   0 1  'uniform'        'sample'  loss []
            2   0 2  'uniform'        'sample'  loss []
            3   1 0  'uniform'        'sample'  loss []
            
            4   0 0  'linear to mean' 'sample'  loss []
            5   0 0  'per condition'  'sample'  loss []
            6   0 0  'per bin'        'sample'  loss []
            
            10  0 0  'uniform'        'shrink'  loss  {exp(-5:0.04:0),exp(-10:0.02:0)}
            11  0 1  'uniform'        'shrink'  loss  {exp(-5:0.04:0),exp(-10:0.02:0)}
            12  0 2  'uniform'        'shrink'  loss  {exp(-5:0.04:0),exp(-10:0.02:0)}
            13  1 0  'uniform'        'shrink'  loss  {exp(-5:0.04:0),exp(-10:0.02:0)}
            
            14  0 0  'linear to mean' 'shrink'  loss  {exp(-5:0.04:0),exp(-10:0.02:0)}
            15  0 0  'per condition'  'shrink'  loss  {exp(-5:0.04:0),exp(-10:0.02:0)}
            16  0 0  'per bin'        'shrink'  loss  {exp(-5:0.04:0),exp(-10:0.02:0)}
            
            })
        end
    end
end