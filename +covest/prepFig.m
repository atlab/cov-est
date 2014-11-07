function prepFig(sz,fname)
f = gcf;
f.PaperUnits = 'centimeters';
f.Units = 'centimeters';
f.PaperSize = sz;
f.PaperPosition = [0 0 sz];
f.Position = [0 0 sz];
a = gca;
a.TickDir = 'out';
a.LabelFontSizeMultiplier = 1.0;
a.FontSize = 8;
a.LineWidth = 0.25;
a.Box = 'off';
warning('off', 'MATLAB:print:CustomResizeFcnInPrint')
s.format = 'eps';
s.Renderer = 'painters';
s.Resolution = 'auto';
%print('-dpng','-r600',fname)
hgexport(f,fname,s)
end
