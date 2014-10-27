%{
sim.Model (lookup) # sampling model 
model :  tinyint   #  model number
-----
model_name : enum('gauss','ising')
%}

classdef Model < dj.Relvar
    methods
        function fill(self)
            self.insert({
                1   'gauss'
                2   'ising'
                })
        end
    end
end