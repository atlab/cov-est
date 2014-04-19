%{
covest.OriTuning (computed) # von Mises tuning for all cells
-> covest.Traces
-----
nshuffles    : smallint  # number of shuffles
von_p_value  : blob      # von mises shuffle p-values
von_r2       : blob      # R-squared of response
von_pref     : blob      # von mises preferred direction
von_base     : blob      # von mises bases
von_amp1     : blob      # von mises preferred amplitude
von_amp2     : blob      # von mises anti-preferred amplitude
von_sharp    : blob      # von mises sharpness
%}

classdef OriTuning < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = covest.Traces & 'ndirs>=8';
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            disp 'loading traces...'
            [X,evokedBins] = fetch1(covest.Traces & key, 'trace_segments','evoked_bins');
            X = permute(squeeze(sum(X(1:min(evokedBins,end),:,:,:))),[3 1 2]);  % discard peristimulus period
            
            disp 'computing von Mises tuning...'            
            nShuffles = 10000;
            
            [von, r2, p] = covest.lib.VonMises2.computeSignificance(X, nShuffles);
            key.nshuffles = nShuffles;
            key.von_r2 = single(r2);
            key.von_base = single(von.w(:,1));
            key.von_amp1 = single(von.w(:,2));
            key.von_amp2 = single(von.w(:,3));
            key.von_sharp= single(von.w(:,4));
            key.von_pref = single(von.w(:,5));
            key.von_p_value = single(p);
            
            self.insert(key)
        end
    end
end
