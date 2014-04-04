function scatter

cove = fetch(pop.AodCovEstimate & 'cov_estim_num=5' & 'bin_opt=0','*');

degrees = arrayfun(@(c) avgNodeDegree(c.sparse), cove);
ranks = arrayfun(@(c) rank(c.lowrank), cove);
scatter(ranks,degrees)

cove = fetch(pop.AodCovEstimate & 'cov_estim_num=2' & 'bin_opt=0','*');
ranks = arrayfun(@(c) rank(c.lowrank), cove);
scatter(ranks,degrees)

end

function degree = avgNodeDegree(C)
p = size(C,1);
degree = sum(sum(~~C.*(1-eye(p))))/p;
end