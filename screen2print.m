%Generate a Bitmap and Vector output from the current figure
%png and eps2 format
%screen2print(filename,OutputSize) [do not specify extension!]
%OutputSize - 300 - HiRes, 100 - normal, 50 -  lowRes
%
% LKB(c) based on Sean P. McCarthy code

function screen2print(filename,OutputSize)

if nargin < 1
error('Not enough input arguments!')
end

oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches','PaperPosition',newpos)

set(gcf, 'renderer', 'painters');
%print bitmap - png
print(gcf, '-dpng', [filename '.png'], ['-r' int2str(OutputSize)]);
%print vector - eps. not recommended if you got a lot of points/lines but useful for any image post-processing
%print(gcf, '-depsc2', [filename '.eps'], ['-r' int2str(OutputSize)]);

drawnow
set(gcf,'Units',oldscreenunits,'PaperUnits',oldpaperunits,'PaperPosition',oldpaperpos)