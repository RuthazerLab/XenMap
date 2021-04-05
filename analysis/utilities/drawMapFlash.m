%{
drawMapFlash.m

input: grid map, lum mask (2D, dimensions must match),color axis limits
(phase),image dimensions
output: image handle
%}

function img = drawMapFlash(imgnum,pixelmap,lummask,flipcmap,clims,imgdim)

img = figure(imgnum);

set(img,'name','Optimal stim position','position',[100 200 512 512],'Color',[1 1 1]);
set(gcf,'units','pixels');
set(gca,'units','pixels');
w_pos = get(gcf, 'position');
set(gca, 'position', [0 0 w_pos(3) w_pos(4)]);

set(img,'Color',[0 0 0]);
img.InvertHardcopy = 'off';

stdadjust = mat2gray(lummask);
Low_High = stretchlim(stdadjust,[0.3,0.9]);
stdadjust = imadjust(stdadjust,Low_High,[]);

posqmap = transpose(reshape(pixelmap,[512,512]));
mapimg = imagesc(posqmap,[1,5]);
hold on;

mapimg.AlphaData = stdadjust;
mapimg.AlphaDataMapping = 'scaled';

hold off;
axis off;

if flipcmap
    colormap(flipud(jet));
else
    colormap(jet);
end
end