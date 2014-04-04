%{
pop.BinOpt (lookup) # ways to obtain binned traces

bin_opt  : tinyint    # binning option
-----
bin_ms   : float   # bin duration
psth_ms  : float   # PSTH duration to subtract from signal before binning
activity_thresh : float  # as fraction of 90th percentile variance
%}

classdef BinOpt < dj.Relvar
    
    properties(Constant)
        table = dj.Table('pop.BinOpt')
    end
    
    methods
        function fill(self)
            s = cell2struct({...
                0     150       1200        0.02
                1     200       1200        0.02
                2     250       1200        0.02
                3     400       1200        0.02                
                },{...
                'bin_opt' 'bin_ms' 'psth_ms' 'activity_thresh'
                },2);
            
            self.inserti(s)
        end
    end
end