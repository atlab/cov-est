%{
pop.TraceOpt (lookup) # trace binning options
trace_opt : tinyint  # trace binning option
-----
bin_width :  float  # (ms) requested bin width
%}

classdef TraceOpt < dj.Relvar
end