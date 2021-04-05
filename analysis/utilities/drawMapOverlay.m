%{
drawMapOverlay.m

input: avg projection,pseudocolor map(overlay),colormap,alpha map(overlay),
color axis limits(overlay),image dimensions
output: image handle

%}
function img = drawMapOverlay(imgnum,stdadjust,img_map,alpha_map,cmap,clims,imgdim)

if isempty(imgdim)
    imgwidth = 512;
    imgheight = 512;
else
    imgwidth = imgdim(1);
    imgheight = imgdim(2);
end


if length(snrmap)>imgwidth
    snrmap = transpose(reshape(snrmap,[imgwidth,imgheight]));
end


img = figure(imgnum);
set(img3,'name','Cell average op stim','position',[150 250 512 512],'Color',[1 1 1]);
set(img3,'Color',[0 0 0]);
img3.InvertHardcopy = 'off';
set(gcf,'units','pixels');
w_pos = get(gcf, 'position');
ax1 = axes;
set(ax1,'units','pixels');
set(ax1, 'position', [0 0 w_pos(3) w_pos(4)]);
ax1_pos = ax1.Position;
meanimg = imagesc(stdadjust);
colormap(ax1,gray);

ax2 = axes;
set(ax2,'Color','none','visible','off');
set(ax2,'units','pixels');
set(ax2, 'position', [0 0 w_pos(3) w_pos(4)]);
cellmapimg = imagesc(transpose(img_map));
cellmapimg.AlphaData = transpose(alpha_map);
cellmapimg.AlphaDataMapping = 'scaled';
if flipcmap
    colormap(ax2,flipud(jet));
else
    colormap(jet);
end

axis off;