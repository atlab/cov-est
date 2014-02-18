function ix = sortcov(C)
% sort the correlation matrix to move highest correlations close to diagonal

T = 1e3;
p = size(C,1);
ix = 1:p;
L = [];
for t = exp(-10*(0:T)/T)
    i = randi(p);
    j = mod(randi(p-1)+i-1,p)+1;
    % evaluate improvement of switching ix(i) and ix(j)
    before = cost(i) + cost(j);
    [ix(i),ix(j)] = deal(ix(j),ix(i));
    improvement = before - (cost(i) + cost(j));
    if improvement + 0.5*t*rand > 0
        L(end+1) = improvement;
    else
        L(end+1) = 0;
        [ix(i),ix(j)] = deal(ix(j),ix(i));
    end
end


    function c = cost(i)
        c = abs((1:p)-i) * (C(ix,ix(i)));
    end
end
