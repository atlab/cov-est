function [avgCorr,avgTunedCorr,avgPCorr,avgTunedPCorr,numTuned,numPos,numNeg,key] = ori(scanNum,oriBounds)
%    oriBounds = [0 10 45 90];

alpha = 0.05;
key = fetch(pop.AodBinnedTraces2 & ...
    sprintf('mod(aod_scan_start_time,1e4)=%d',scanNum) & 'bin_opt=0');
key.cov_estim_num=14;

cove = fetch(pop.AodCovEstimate2*pro(pop.AodBinnedTraces2,'cellnums') & key, '*');
sz = size(cove.sparse);
[ix,jx] = ndgrid(1:sz(1),1:sz(2));
offDiag = jx>ix;
offDiag = offDiag(:);
oris = fetch(pop.AodVonMises & key, '*');
Sigma0 = fetch1(pop.AodCovEstimate2 & setfield(key,'cov_estim_num',0),'cov_matrix'); %#ok<SFLD>
C0 = corrcov(Sigma0);
avgCorr = mean(C0(offDiag));
P0 = -corrcov(cove.sparse-cove.lowrank*cove.lowrank');
avgPCorr = mean(P0(offDiag));

if isempty(oris)
    s = [1,length(oriBounds)-1];
    avgTunedCorr = nan(s);
    avgTunedPCorr = nan(s);
    numTuned = zeros(s);
    numPos   = zeros(s);
    numNeg   = zeros(s);
else
    oris.von_pref = mod(oris.von_pref*180/pi,360);
    deltaOri = nan(size(cove.sparse));
    p = size(cove.sparse,1);
    for i=2:p
        ix = find(oris.cellnums==cove.cellnums(i));
        assert(length(ix)==1)
        for j=1:i-1
            jx = find(oris.cellnums==cove.cellnums(j));
            assert(length(jx)==1)
            if oris.von_p_value(ix) < alpha && oris.von_p_value(jx) < alpha
                d = ne7.rf.oriDiff(oris.von_pref(ix),oris.von_pref(jx));
                deltaOri(i,j) = d;
                deltaOri(j,i) = d;
            end
        end
    end
    tunedPairs = offDiag & ~isnan(deltaOri(:));
    avgTunedCorr = arrayfun(@(i) mean(C0(offDiag & between(deltaOri(:), oriBounds(i+(-1:0))))), 2:length(oriBounds));
    avgTunedPCorr = arrayfun(@(i) mean(P0(offDiag & between(deltaOri(:), oriBounds(i+(-1:0))))), 2:length(oriBounds));
    numTuned = arrayfun(@(i) sum(tunedPairs & between(deltaOri(:), oriBounds(i+(-1:0)))), 2:length(oriBounds));
    numPos   = arrayfun(@(i) sum(cove.sparse(tunedPairs & between(deltaOri(:), oriBounds(i+(-1:0)))) <0), 2:length(oriBounds));
    numNeg   = arrayfun(@(i) sum(cove.sparse(tunedPairs & between(deltaOri(:), oriBounds(i+(-1:0)))) >0), 2:length(oriBounds));
end
end

function yes = between(a,rng)
yes = a > rng(1) & a<=rng(2);
end