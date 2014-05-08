classdef VonMises2 < handle
    % VonMises2 - fast fitting of vonMises tuning curves.
    % For usage example, execute and examine covest.lib.VonMises2.test
    %
    % Dimitri Yatsenko, 2012-07-09
    
    properties(SetAccess = private)
        w  % nx5 matrix of coefs in equation
        %  w1 + w2 * exp(c-1) + w3 * exp(-c-1)    with  c = w4 * cos(phi - w5)
        %  where  w1 >= w2 >= 0
        phi
        nDirs
    end
    
    methods(Static)
        function test(n)
            % synthesize n randomized tuning curves
            if ~nargin, n = 5;  end
            
            % generate random tuning functions
            ph = (0:22.5:359)/180*pi;
            base = (1:n)'*2;
            smallPeak = random('gam',1,1,n,1);
            bigPeak = smallPeak + random('gam',1,1,n,1);
            sharpness = random('gam',1,4,n,1)+1;
            prefDir = random('uni',0,2*pi,n,1);
            v = [base bigPeak smallPeak sharpness, prefDir];
            x = compute(covest.lib.VonMises2(v), ph);
            x = x+0.2*random('norm', 0, 1, size(x));  % add noise
            
            % now fit the tuning functions
            tic
            f = fit(covest.lib.VonMises2, x');
            toc
            
            if true
                % if eight traces or fewer, print and plot results
                disp 'Original coeffs:'
                disp(v)
                plot(ph*180/pi, x, 'o', 'MarkerSize', 10)
                hold on
                ph = (0:5:360)/180*pi;
                plot(ph*180/pi, f.compute(ph));
                xlabel 'direction (degrees)'
                disp 'Fitted coeffs:'
                disp(f.w)
                xlim([0 360])
                grid on
                hold off
            end
        end
        
        function [von, r2, p] = computeSignificance(responses, nShuffles)
            % compute von mises tuning and its siginifiance by shuffling
            % responses must be a matrix of responses with dimensions
            %   nCells * nOris * nTrials
            % and it may contain NaNs
            [von,r2] = computeTuning(responses);
            
            sz = size(responses);
            responses = reshape(responses, sz(1),[]);
            p = 0.5/nShuffles;
            for i=1:nShuffles
                if ~mod(i,250) || any(i==[1 nShuffles])
                    fprintf('Shuffles [%4d/%4d]\n', i, nShuffles)
                end
                [~, r2_] = computeTuning(reshape(responses(:,randperm(end)),sz));
                p = p + (r2_>=r2)/nShuffles;
            end
            
            function [von, r2] = computeTuning(x)
                von = fit(covest.lib.VonMises2, nanmean(x,3)');
                y = bsxfun(@minus, x, von.compute);
                r2 = 1-nanvar(reshape(y,size(x,1),[])')./nanvar(reshape(x,size(x,1),[])');
            end
        end
    end


    
    methods
        function self = VonMises2(w)
            if nargin
                assert(isnumeric(w) && size(w,2)==5, 'invalid input')
                self.w = w;
            end
        end

        
        function self = fit(self, F)
            % fits data F = k*n matrix where n is the number of tuning curves
            % and k is the number of uniformly distributed directions
            
            % subtract the mean
            F = double(F);
            self.nDirs = size(F,1);
            assert(~bitand(self.nDirs,1), 'Must have an even number of sampled angles')
            M = mean(F);
            F = bsxfun(@minus, F, M);
            
            % find initial solutions exhaustively among several preferred
            % directions at mid-level sharpness
            self.phi = (0:self.nDirs-1)/self.nDirs*2*pi;
            thetaStep = pi/7;
            sharp = 5;
            for theta = 0:thetaStep:(pi-1e-5)
                [ww,r2] = self.fitLinear(F, ...
                    repmat([0 nan nan sharp theta], size(F,2), 1));
                if ~theta
                    bestR2 = r2;
                    self.w = ww;
                else
                    ix = r2 > bestR2;
                    self.w(ix,:) = ww(ix,:);
                    bestR2(ix) = r2(ix);
                end
            end
            
            thetaStep = pi/12;
            sharpStep = log(4);  % in natural log space
            minSharp = 0.7;
            
            for iter=1:3
                % adjust theta.  2-3 iterations are sufficient
                ww = self.w;
                ww(:,5) = self.w(:,5) - thetaStep;
                [ww,left] = self.fitLinear(F,ww);
                ww(:,5) = self.w(:,5) + thetaStep;
                [ww,right] = self.fitLinear(F,ww);
                adjust = (right-left)./(2*bestR2-left-right)/2;
                ww(:,5) = self.w(:,5) + thetaStep*max(-1, min(1, adjust));
                [ww,r2] = self.fitLinear(F,ww);
                ix = r2>=bestR2;
                self.w(ix,:) = ww(ix,:);
                bestR2(ix) = r2(ix);
                thetaStep = thetaStep/10;
                
                % adjust sharpness (binary search, need more iterations)
                for j=1:1+iter
                    for sgn=[-1 1]
                        ww = self.w;
                        ww(:,4) = (self.w(:,4)-minSharp).*exp(sgn*sharpStep)+minSharp;
                        [ww,r2] = self.fitLinear(F,ww);
                        ix = r2>=bestR2;
                        self.w(ix,:) = ww(ix,:);
                        bestR2(ix) = r2(ix);
                        sharpStep = sharpStep*0.75;
                    end
                end
            end
            
            % flip by 180 degrees if peak2 > peak1
            ix = self.w(:,2) < self.w(:,3);
            self.w(ix,2:3) = self.w(ix,[3 2]);
            self.w(ix,5) = self.w(ix,5) + pi;
            self.w(:,5) = mod(self.w(:,5),2*pi);
            
            % restore base
            mm = mean(self.makePeak(self.w),2);
            self.w(:,1) = M' - (self.w(:,2)+self.w(:,3)).*mm;
        end
        
        function F = compute(self, phi)
            if nargin < 2
                phi = self.phi;
            end
            g2 = cos(bsxfun(@minus, phi, self.w(:,5)));
            g1 = exp(bsxfun(@times, self.w(:,4),  g2-1));
            g2 = exp(bsxfun(@times, self.w(:,4), -g2-1));
            F = bsxfun(@times, self.w(:,2), g1) + ...
                bsxfun(@times, self.w(:,3), g2);
            F = bsxfun(@plus, self.w(:,1), F);
        end
        
               
    end
    

    
    methods(Access = private)
        
        function [ww, r2] = fitLinear(self, F, ww)
            % update the linear coefficients in ww
            [g1,g2] = self.getPeaks(ww);
            
            % covariance matrix [a c; c b]
            a = sum(g1.*g1,2);
            b = sum(g2.*g2,2);
            c = sum(g2.*g1,2);
            
            % invert covariance matrix
            idet = 1./(a.*b - c.*c);
            b = idet.*b;
            c = idet.*c;
            
            % solve for linear coefficients
            h1 =  bsxfun(@times, b, g1) - bsxfun(@times, c, g2);
            h2 =  bsxfun(@times, b, g2) - bsxfun(@times, c, g1);
            ww(:,2) = max(0,sum(h1'.*F)');
            ww(:,3) = max(0,sum(h2'.*F)');
            
            % compute sum of squares to measure fit quality
            r2 = sum((...
                bsxfun(@times, ww(:,2), g1) + ...
                bsxfun(@times, ww(:,3), g2)).^2,2);
        end
        
        
        function g = makePeak(self, w)
            g = exp(bsxfun(@times, w(:,4), ...
                cos(bsxfun(@minus, self.phi, w(:,5)))-1));
        end
        
        
        function [g1, g2] = getPeaks(self, w)
            g1 = self.makePeak(w);
            g1 = bsxfun(@minus, g1, mean(g1,2));  % make zero-mean
            g2 = circshift(g1,[0 self.nDirs/2]);
        end
        
    end
end