%{
pop.CovEstimator (lookup) # covariance table estimators
cov_estim_num   :  tinyint   #  covariance estimator number
-----
estimator_name              : varchar(80)                   # name for legends, etc.
estimator_call              : varchar(80)                   # identifier for switch cases
hyperparam_space            : longblob                      # a cell array with arrays of values defining the search space
loss_fun                    : varchar(255)                  # inline function
%}

classdef CovEstimator < dj.Relvar
    
    properties(Constant)
        table = dj.Table('pop.CovEstimator')
    end
    
    methods
        function fill(self)
            loss =  '@(S,Sigma)(trace(Sigma/S)+logDet(S))/size(S,1)';
            
            searchSpace = {};
            self.inserti(struct('cov_estim_num',0,'estimator_name','sample','estimator_call','sample','hyperparam_space',{searchSpace},'loss_fun',loss))
            
            searchSpace = {1,exp(-6:0.02:0)};
            self.inserti(struct('cov_estim_num',11,'estimator_name','shrink','estimator_call','shrink','hyperparam_space',{searchSpace},'loss_fun',loss))
            
            searchSpace = {0:40,exp(-6:0.02:0)};
            self.inserti(struct('cov_estim_num',12,'estimator_name','multifactor','estimator_call','factor','hyperparam_space',{searchSpace},'loss_fun',loss))
            
            searchSpace = {0,exp(-6:0.02:0)};
            self.inserti(struct('cov_estim_num',13,'estimator_name','sparse+lowrank','estimator_call','lv-glasso','hyperparam_space',{searchSpace},'loss_fun',loss))
            
            searchSpace = {0:40,exp(-6:0.02:0)};
            self.inserti(struct('cov_estim_num',14,'estimator_name','sparse+lowrank','estimator_call','lv-glasso','hyperparam_space',{searchSpace},'loss_fun',loss))

            searchSpace = {exp(-4:0.04:0),exp(-6:0.02:0)};
            self.inserti(struct('cov_estim_num',15,'estimator_name','sparse+lowrank','estimator_call','lv-glasso','hyperparam_space',{searchSpace},'loss_fun',loss))

            searchSpace = {exp(-5:0.04:0),exp(-10:0.02:0)};
            self.inserti(struct('cov_estim_num',21,'estimator_name','shrink','estimator_call','shrink','hyperparam_space',{searchSpace},'loss_fun',loss))

        end
    end
end