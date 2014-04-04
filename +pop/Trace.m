%{
pop.Trace (imported) # my newest table
-> pop.Segment
trace_idx      : smallint  # cell index within the segmentation
---
gtrace         : longblob  # fluorescence trace in green channel
centroid_x     : float  # (pixels) centroid x-coordinate
centroid_y     : float  # (pixels) centroid y-coordinate
mask_pixels    : mediumblob                    # indices of segment pixels in the image
g_brightness   : float        # mean brightness inside cell
g_background   : float        # 40th percentile brightness in surrounding annulus
r_brightness   : float        # mean brightness inside cell
r_background   : float        # 40th percentile brightness in surrounding annulus
%}

classdef Trace < dj.Relvar & dj.AutoPopulate
    
    properties(Constant)
        table = dj.Table('pop.Trace')
        popRel = pop.Segment
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            % compute pixel pitch (um/pixel)
            [px,py] = fetch1(tp.Align & key, ...
                'um_width/px_width->px', 'um_height/px_height->py');
            assert(max(px/py, py/px)<1.1, ...
                'the present algorithm cannot accept non-isometric pixels')
            pitch = (px+py)/2;  %microns per pixel

            mask = logical(fetch1(pop.Segment & key, 'manual_mask'));

            % get mean frames and movie
            [greenImage, redImage] = fetch1(tp.FineAlign & key, 'fine_green_img', 'fine_red_img');
            movie = tp.utils.Movie(key);
            
            % extract traces
            disp 'loading movie...'
            X = movie.getFrames(1,1:movie.nFrames);
            sz = size(X);
            regions = regionprops(bwconncomp(mask,4), 'Centroid','PixelIdxList');
            X = reshape(X,[],sz(3))';
            
            radius = 15/pitch; % 15-micron radius
            
            for iTrace = 1:numel(regions)
                tuple = key;
                tuple.trace_idx = iTrace;
                tuple.gtrace = mean(X(:,regions(iTrace).PixelIdxList),2);
                tuple.centroid_x = regions(iTrace).Centroid(1);
                tuple.centroid_y = regions(iTrace).Centroid(2);
                tuple.mask_pixels = regions(iTrace).PixelIdxList;
                tuple.g_brightness = mean(greenImage(regions(iTrace).PixelIdxList));
                tuple.g_background = getBackground(greenImage, regions(iTrace).PixelIdxList, radius);
                tuple.r_brightness = mean(redImage(regions(iTrace).PixelIdxList));
                tuple.r_background = getBackground(redImage, regions(iTrace).PixelIdxList, radius);
                
                self.insert(tuple)
            end
        end
    end
end



function background = getBackground(img, maskPixels, radius)
% Michelson contrast of pixels at maskPixels compared to others within radius
sz = size(img);
[y,x] = ind2sub(sz,maskPixels);
y = mean(y);
x = mean(x);
[yi,xi] = ndgrid(...
    max(1,round(y-radius)):min(sz(1),round(y+radius)),...
    max(1,round(x-radius)):min(sz(2),round(x+radius)));
hood = (yi-y).^2+(xi-x).^2 < radius^2;
hood = setdiff(sub2ind(sz, yi(hood), xi(hood)), maskPixels);
background = quantile(img(hood),0.4);
end