%{

an_gradaxis.m

Displays 3D phase map and calculates map gradient vector

%}

%% user parameters

plotstyle = 1;      %1 = surface, 2 = scatter
nsamples = 5000;    %for scatter plot: number of points to show
markersize = 2;     %for scatter plot: size of markers
plotline = 1;       %plot gradient axis line
linecolor = [0,0,0];
linetype = '-'; 

imgxflip = 1;       %lr flip images

%% ---end of user parameters---

addpath('exampledata');
addpath('examplemasks');

date = '20200302';
usemask = 1;
map1list = [2036 2120];
map2list = [2041 2125];
goslice = 2:8;
stepsize = 10;
imgwidth=256;

linecolor = linecolor./256;
imgheight = imgwidth;
filt=1;
snrthres = 1;
takesamples = 1;
gaussfilt = 0;

fieldwidth = 254.13;
pixelrat = imgwidth/fieldwidth;
nslices = length(map1list)*length(goslice);
imgdepth =  round(stepsize*(nslices+1));

axislims.x = [0 imgwidth];
axislims.y = axislims.x;
axislims.z = [0 (nslices+1)*stepsize * pixelrat];

scatterXcache = {nslices};
scatterYcache = {nslices};
scatterZcache = {nslices};
scatterCcache = {nslices};
readslice = 1;
for j=1:length(map1list)
    map1id = num2str(map1list(j));
    if map1list(j)<1000
        map1id = ['0' map1id];
    end
    
    map2id = num2str(map2list(j));
    if map2list(j)<1000
        map2id = ['0' map2id];
    end
    
    ID = [map1id '-' map2id];
    
    if j==1 %read this only once
        datafile = ['StimulusData_' date '_' map1id(1:2) '.' map1id(3:4) '.mat'];
        load(datafile);
        pblank = setup.DriftSettings.Blanktime/(setup.DriftSettings.Blanktime+1/setup.DriftSettings.Frequency);
        orientation = setup.DriftSettings.Orientation;
        if strcmp(orientation,'horizontal')
            clims = [-pi -pi+2*pi*(1-pblank)*0.5];
            stimaxis = 'elevation';
            ylims = [-pi,-pi+2*pi*(1-pblank)*0.5];
        else
            clims = [-pi -pi+2*pi*(1-pblank)];
            stimaxis = 'azimuth';
            ylims = [-pi,-pi+2*pi*(1-pblank)];
        end    
    end
    
    
    for crtslice = goslice
        savenotes = ['_sl' num2str(crtslice) '_filt=' num2str(filt)];
        mapfile = [date '_' map1id '-' map2id savenotes '.mat'];
        load(mapfile);
        
        
       %% apply neuropil mask and snr threshold
        volumedat(:,:,readslice) = map_sub;
        crtslicemask = ones(imgwidth,imgheight);
        snrmap = transpose(reshape(snrmap,[imgwidth,imgheight]));
        crtslicemap = volumedat(:,:,readslice);
        
        if usemask==1
            maskfile = ['AVG_Exp_' map2id '_slice' num2str(crtslice) '.tif'];
            inmask = imread(maskfile);
        elseif usemask==2
            maskfile = ['Map_Exp_' map1id '-' map2id '_slice' num2str(crtslice) '.tif'];
            inmask = imread(maskfile);
            inmask_resize = imresize(inmask, [256 256]);
            inmask_bw = im2bw(inmask_resize,0.01);
            inmask = inmask_bw;
        end
        
        if gaussfilt
            crtslicemap = imgaussfilt(crtslicemap,gaussfilt);
        end
        
        crtslicemap(inmask<=0) = NaN;
        crtslicemap(snrmap<snrthres) = NaN;
        crtslicemask(inmask<=0) = 0;
        maskdat2(:,:,readslice) = crtslicemask;
        crtslicemask(snrmap<snrthres) = 0;
        maskdat(:,:,readslice) = crtslicemask;
        volumedat(:,:,readslice) = crtslicemap;
                
        
        %% for scatterplot
        if plotstyle == 2
            scattermesh = 1:imgwidth*imgheight;
            scatterpoints = scattermesh(transpose(crtslicemask)>0);
            if length(scatterpoints)>nsamples
                scatterinds = randsample(scatterpoints,nsamples);
            else
                scatterinds = scatterpoints;
            end

            
            [scatterX,scatterY] = ind2sub([imgwidth,imgheight],scatterinds);
            scatterZ = (ones(1,length(scatterinds))*crtslice + length(goslice)*(j-1)) * stepsize * pixelrat;
            
            map_trans = transpose(map_sub);
            scatterC = crtslicemap(scatterinds);
            if nslices == 1
                writeslice = 1;
            else
                writeslice = crtslice + length(goslice)*(j-1);
            end
            scatterXcache{writeslice}=scatterX;
            scatterYcache{writeslice}=scatterY;
            scatterZcache{writeslice}=scatterZ;
            scatterCcache{writeslice}=scatterC;
        end      
                
        readslice = readslice+1;
    end
    
end


%% now calculate gradient

[dx,dy,dz] = gradient(volumedat);
sumx = nansum(nansum(nansum(dx)));
sumy = nansum(nansum(nansum(dy)));
sumz = nansum(nansum(nansum(dz)));
%z direction need to be corrected for scale difference: neighboring pixels
%in z direction are stepsize*pixelrat pixels apart
sumz_corr = sumz/(stepsize*pixelrat);

vector = [sumx,sumy,sumz_corr];
vector_norm = vector/norm(vector);
npoints = sum(sum(sum(~isnan(volumedat))));
gradsize = norm(vector)/npoints;


%% plotting

%find a centroid point (on a center slice) to pass axis line
centerslice = floor(nslices/2)+1;
centerseg = maskdat2(:,:,centerslice);
centerseg(centerseg>0)=1;
rps = regionprops(centerseg,'centroid');
tcentroid_xy = rps.Centroid;
tcentroid_z = centerslice * stepsize * pixelrat;
tcentroid = [tcentroid_xy(1) tcentroid_xy(2) tcentroid_z];

%calculate line segment for plotting
[lsx,lsy,lsz]=plotvector(vector,tcentroid,axislims);


x=(1:1:imgwidth);
y=(1:1:imgwidth);
z=(1:1:nslices);
[X,Y,Z] = meshgrid(x,y,z);
[XD,YD]=meshgrid(x,y);


img = figure(1);
set(img,'position',[100 100 500 400]);
clf;
hold on;

if plotstyle == 1  %surf plots
    
    for countslice = 1:nslices
        ZD=ones(imgwidth) * countslice * stepsize * pixelrat;
        mapsurf=surface(XD,YD,ZD,volumedat(:,:,countslice),'FaceColor','texturemap','EdgeColor','none','CDataMapping','scaled');
        mapsurf.FaceAlpha='texturemap';
        
        if usemask
            alphamask = maskdat(:,:,countslice);
            alphamask(alphamask~=0)=1;
            mapsurf.AlphaData=alphamask;
            mapsurf.AlphaDataMapping = 'scaled';
        end
        
        alim([0,1]);
    end
    
else %scatter plots
    
    [X,Y,Z] = meshgrid(imgwidth,imgheight,imgdepth);
    scatter3(scatterXcache{1},scatterYcache{1},scatterZcache{1},markersize,scatterCcache{1},'filled');
    hold on;
    for crtslice = 2:nslices
        scatter3(scatterXcache{crtslice},scatterYcache{crtslice},scatterZcache{crtslice},markersize,scatterCcache{crtslice},'filled');
    end   
    xlim([0,imgwidth]);
    ylim([0,imgheight]);
    zlim([0,imgdepth]);
    set(gca,'DataAspectRatio',[1 1 1])
end



set(gca,'Zdir','reverse');
set(gca, 'YDir','reverse');
if imgxflip
    set(gca, 'XDir','reverse');
end

colormap(jet);
caxis(clims);


if plotline == 1
	plot3(lsx,lsy,lsz,linetype,'color',linecolor,'linewidth',2);
end

xlim(axislims.x);
ylim(axislims.y);
zlim(axislims.z);

xlim([1,256]);
ylim([1,256]);

zlim(axislims.z);


set(gca,'box','on','linewidth',1.5);

if imgxflip
    view(-11,10)
else
    view(78,14)
end

grid off;
hold off;

gradvector = [sumx sumy sumz_corr];
%set(gca,'DataAspectRatio',[1 1 1])

disp(['Gradient vector: (' num2str(sumx) ' ' num2str(sumy) ' ' num2str(sumz_corr) ')'])
