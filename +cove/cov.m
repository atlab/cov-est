% covariance matrix without mean subtraction

function [C,n] = cov(X)
ix = ~any(isnan(X),2);
n = sum(ix);
C = X(ix,:)'*X(ix,:)/(n-1);