function g = doppler(m)
% doppler(m)
%   doppler(M) returns an M-by-3 matrix containing a doppler colormap.
%   doppler, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%             colormap(doppler)
%
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG, COLORMAP, RGBPLOT.
%
% Dimitri Yatsenko, 2011-08-10

if nargin < 1, m = size(get(gcf,'colormap'),1); end
m  = floor(m/2)-1;
ix = [-m:-1 0 0 1:m]'/m;
g  = [1-max(0,-ix) 1-max(ix,-ix) 1-max(0,ix)];

%z = exp(-abs(-1:1/31.5:1)'.^1.0);
%g = bsxfun(@times, g.^1.0, z);