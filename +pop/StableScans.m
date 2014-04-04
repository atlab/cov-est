%{
pop.StableScans (manual) # high-quality scans selected for processing

-> acq.AodScan
---
%}


classdef StableScans < dj.Relvar
    
    properties(Constant)
        table = dj.Table('pop.StableScans')
    end
end