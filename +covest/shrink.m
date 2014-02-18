function [C,lambda] = shrink(X,lambdaVar,lambda)
% shrink variance toward median by amount lambdaVar [0,1]
% shrink covariance toward diagonal by amount lambda [0,1]
%
% If lambdaVar or lambda are not provided, they are estimated analytically as
% implemented by the function cov.shrink of the corpcor R library.
%
% See also Schafer & Strimmer (2005) "A Shrinkage Approach to Large-Scale
% Covariance Matrix Estimation and Implications for Functional Genomics"
% http://uni-leipzig.de/~strimmer/lab/publications/journals/shrinkcov2005.pdf
%
% INPUTS:
%    X - n*p matrix with n p-dimensional observations. Assumed zero mean.
%    lambdaVar (optional) - variance shrinkage intensity toward median
%    lambda (optional) - covariance shrinkage toward diagonal
%
% OUTPUTS:
%    C - covariance estimate
%    lambda - covariance shrinkage intensity

% Dimitri Yatsenko, July 2013


[n,p] = size(X);

% shrink variance toward median variance
vars = mean(X.^2);
medVar = median(vars);

doLambdaVar  = nargin < 2 || isempty(lambdaVar) || isnan(lambdaVar);
doLambda = nargin < 3;


if doLambdaVar || doLambda
    % estimate
    SS = nan(n,p,p);
    for i=1:n
        SS(i,:,:) = X(i,:)'*X(i,:);
    end
    V = squeeze(var(SS))/n;
    M = squeeze(mean(SS));
end

if doLambdaVar    
    lambdaVar = min(1,sum(diag(V))./sum(diag(M)-medVar).^2);
end

target = diag(lambdaVar*medVar + (1-lambdaVar)*vars);
if doLambda
    % estimate optimal covariance shrinkage
    lambda = min(1, (sum(V(:))-lambdaVar*sum(diag(V)))./sum(sum((M-target).^2)));
end

% shrink sample covariance toward diagonal
C = (1-lambda)*(X'*X/n) + lambda*target;