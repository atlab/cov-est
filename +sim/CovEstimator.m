%{
sim.CovEstimator (lookup) # my newest table
cov_estim_num   :  tinyint   #  covariance estimator number
-----
estimator_name              : varchar(80)                   # name for legends, etc.
estimator_call              : varchar(80)                   # identifier for switch cases
hyperparam_space            : longblob                      # a cell array with arrays of

%}

classdef CovEstimator < dj.Relvar
    
    properties(Constant)
        table = dj.Table('sim.CovEstimator')
    end
    
    methods
        function fill(self)
            
            searchSpace = {};
            self.inserti(struct('cov_estim_num',0,'estimator_name','sample','estimator_call','sample','hyperparam_space',{searchSpace}))
            
            searchSpace = {exp(-5:0.04:0),exp(-10:0.02:0)};
            self.inserti(struct('cov_estim_num',1,'estimator_name','shrink','estimator_call','shrink','hyperparam_space',{searchSpace}))
            
            searchSpace = {1:20,exp(-10:0.02:0)};
            self.inserti(struct('cov_estim_num',2,'estimator_name','multifactor','estimator_call','factor','hyperparam_space',{searchSpace}))
            
            searchSpace = {0,exp(-10:0.02:0)};
            self.inserti(struct('cov_estim_num',3,'estimator_name','sparse','estimator_call','lv-glasso','hyperparam_space',{searchSpace}))
            
            searchSpace = {0:20,exp(-10:0.02:0)};
            self.inserti(struct('cov_estim_num',4,'estimator_name','sparse+lowrank','estimator_call','lv-glasso','hyperparam_space',{searchSpace}))

            searchSpace = {0:20,exp(-10:0.02:0)};
            self.inserti(struct('cov_estim_num',5,'estimator_name','sparse+lowrank','estimator_call','lv-glasso','hyperparam_space',{searchSpace}))
            
        end
    end
end
