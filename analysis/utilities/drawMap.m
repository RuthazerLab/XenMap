%{
drawMap.m
read phase map data and display phase map (masked by SNR)

input: phase map, snr mask (2D, dimensions must match),color axis limits
(phase),image dimensions
output: image handle

%}
function img = drawMap(imgnum, map_sub,snrmap,clims,imgdim)

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
set(img,'name','Phase map','position',[100 200 512 512],'Color',[1 1 1]);
set(gcf,'units','pixels');
set(gca,'units','pixels');
w_pos = get(gcf, 'position');
set(gca, 'position', [0 0 w_pos(3) w_pos(4)]);
img.InvertHardcopy = 'off';
mapimg = imagesc(map_sub,clims);
axis off;
set(img,'Color',[0 0 0]);
snrmask = mat2gray(snrmap);
Low_High = stretchlim(snrmask,[0.3,0.99]);
snrmask = imadjust(snrmask,Low_High,[]);
hold on;
mapimg.AlphaData = snrmask;
mapimg.AlphaDataMapping = 'scaled';
hold off;
colormap jet;
end