% covariance matrix without mean subtraction

function C = cov(X)
C = X'*X/size(X,1);