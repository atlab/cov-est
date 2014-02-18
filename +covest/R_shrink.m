function C = R_shrink(X)
% use R's corpcor package to do Ledoit-Wolf shrinkage

% write out the algorithm output
inFile = '~/glasso_input.dat';
outFile = '~/glasso_output.dat';
f = fopen(inFile,'w');
p = size(X,2);
fwrite(f,[p X(:)'],'double');
fclose(f);

% run the R glasso routine
Rdir = [fileparts(mfilename('fullpath')) '/../../R'];
command = sprintf('exec -c Rscript %s/runShrink.R %s %s',Rdir, inFile, outFile);
system(command);

% read the result
f = fopen(outFile,'r');
lambdas = fread(f,'double');
fclose(f);

% apply shrinkage to variances and correlations
lambda = lambdas(1);
lambdaVar = lambdas(2);
C = covest.cov(X);
sigma = sqrt(diag(C));
R = bsxfun(@rdivide,bsxfun(@rdivide,C,sigma),sigma');
R = (1-lambda)*R + lambda*eye(p);
sigma2 = sigma.^2;
sigma2 = (1-lambdaVar)*sigma2 + lambdaVar*median(sigma2);
sigma = sqrt(sigma2);
C = bsxfun(@times,bsxfun(@times,R,sigma),sigma');