%{
pop.Segment (imported) # manual segmentaion using different kinds of source images
-> tp.FineAlign
-> pop.SegOpt
---
manual_mask  :  longblob   # binary mask on the finely aligned image
%}

classdef Segment < dj.Relvar & dj.AutoPopulate
    
    properties(Constant)
        table = dj.Table('pop.Segment')
        popRel = tp.FineAlign * pop.SegOpt;
    end
    
    methods(Access=protected)
          function makeTuples(self, key)
            bw = pop.Segment.outlineCells(key,[]);
            
            assert(~isempty(bw), 'user aborted segmentation')
            key.manual_mask = bw;
            self.insert(key)
        end
        
        
        function redo(self)
            % edit existing segmentations
            for key = fetch(self)
                bw = fetch1(tp.SegmentManual & key, 'manual_mask');
                bw = tp.SegmentManual.outlineCells(key, bw);
                if ~isempty(bw)
                    del(tp.Segment & key)
                    insert(tp.Segment, key)
                    key.manual_mask = bw;
                    insert(tp.SegmentManual, key)
                    disp 'updated mask'
                end
            end
        end
    end
    
    
    methods(Static,Access=private)
        
        function bw = outlineCells(key,bw)
            opt = fetch(pop.SegOpt & key, '*');
            switch opt.source_image
                case 'green'
                    img = fetch1(tp.FineAlign & key, 'fine_green_img');
                    img = sqrt(img);
                    img = max(0,img-quantile(img(:),0.01));
                    img = min(1,img/quantile(img(:),0.995));
                    
                case 'fine orimap'
                    cond = struct('ca_opt', 7);
                    img = fetch1(tp.FineVonMap & key & cond, 'von_r2');
                    img = sqrt(img);
                    mx = quantile(img(:),0.995);
                    fprintf('Max = %f\n', mx);
                    mx = 0.5;   % always scale the same way
                    img = 1-max(0, min(1, img/mx));
            end
            f = figure;
            imshow(img)
            set(gca, 'Position', [0.05 0.05 0.9 0.9]);
            pos = get(f, 'Position');
            set(f, 'Position', [pos(1:2)/4 pos(3:4)*4])
            bw = ne7.ui.drawCells(bw);
            close(f)
        end
    end
end
