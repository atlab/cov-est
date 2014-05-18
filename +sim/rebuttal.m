classdef rebuttal
    
    methods(Static)
        function droppedUnit
            
            % produce the true covariance matrix C with sparse interactions
            p = 50;
            n = 5*p;
            rng(1)
            V =  lognrnd(1,.25,p,1);  % variances
            X = mvnrnd(zeros(1,p),diag(V),n);
            trueC = cove.cov(X);
            cove.set('refit',true)
            res = cove.lvglasso(trueC,0.2,1.0,cove.set);
            
            % remove the 5 most connected unit
            nOmit = 5;
            d = sum(logical(res.S));
            [~,ix] = sort(-d);
            ix = ix(nOmit+1:end);
            trueC = inv(res.S(ix,ix));
            altC = inv(res.S);
            altC = altC(ix,ix);
            p = p - nOmit;
            
            % draw a large sample from each distribution
            n = 4000;
            
            rng(2)
            X = mvnrnd(zeros(1,p),trueC,n);  % truth
            X = reshape(X,1,1,size(X,1),size(X,2));
            rng(2)
            Y = mvnrnd(zeros(1,p),altC,n); % dropped unit
            Y = reshape(Y,1,1,size(Y,1),size(Y,2));
            
            cove.set('refit',true)
            
            % infer models by cross-validation
            loss = @(S,Sigma)(trace(Sigma/S)+cove.logDet(S))/size(S,1);
            hyperSpace =  {0,exp(-6:.05:0),exp(-6:.05:1)};
            reg = 'lv-glasso';
            [hypers1,visited1,losses1] = cove.crossEstimateHyper(X, 1, loss, reg, hyperSpace);
            [C1,~,extras1] = cove.estimate(X, [], 1, reg, hypers1);
            
            [hypers2,visited2,losses2] = cove.crossEstimateHyper(Y, 1, loss, reg, hyperSpace);
            [C2,~,extras2] = cove.estimate(Y, [], 1, reg, hypers2);
            
            disp done
            
        end
    end
    
end