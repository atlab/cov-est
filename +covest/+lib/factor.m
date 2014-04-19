% function [L,Ph,LL]=factor(X,K,cyc,tol);
%
% Fast Maximum Likelihood Factor Analysis using EM
%
% X - data or covariance (if square) matrix
% K - number of factors
% cyc - maximum number of cycles of EM (default 100)
% tol - termination tolerance (prop change in likelihood) (default 0.0001)
%
% L - factor loadings
% Ph - diagonal uniquenesses matrix
% LL - log likelihood curve
%
% Iterates until a proportional change < tol in the log likelihood
% or cyc steps of EM.
%
% Adopted by Dimitri Yatsenko from Zoubin Gharamani


function [L,Ph,LL]=factor(X,K,cyc,tol)

if nargin<4, tol=1e-9; end;
if nargin<3, cyc=200; end;

[n,p] = size(X);   assert(n>=p)
if n>p, X = cov(X);  end   % if not square, compute covariance

[L,d] = eigs(X,K); L = L*sqrt(d);  % initialize with PCA
Ph=diag(X);  I=eye(K);  LL=nan(1,cyc);

for i = 1:cyc
    % Estimation
    LP=bsxfun(@rdivide,L,Ph);
    M=diag(1./Ph)-LP/(I+L'*LP)*LP';  % inverse covariance
    b=L'*M;  Xb=X*b';
    
    % Maximization
    L=Xb/(I+b*(Xb-L));
    Ph=diag(X)-diag(L*Xb');
    
    % compute log likelihood and check termination conditions
    [C,p] = chol(M);
    if p
        LL(i) = inf;
    else
        LL(i) = n/2*(-p*log(2*pi)+2*sum(log(diag(C)))-sum(sum(M.*X)));
    end
    if isinf(LL(i)) || i>2 && LL(i)-LL(i-1)<tol*(LL(i)-LL(1)), break, end
end