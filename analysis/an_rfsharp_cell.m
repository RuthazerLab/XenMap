%{

an_rfsharp_cell.m

Calculates receptive field sharpness for cell bodies
(RF sharpness = average ?F/F? response to the 2 or 3 stimulus positions
closest to the optimal stimulus position (grid) divided by the average 
response to the remaining stimulus positions in the periphery)

%}


%% user parameters

cellrthres = 2; %threshold for peak response amplitude (exclude cells with peak response below threshold)
imgxflip = 1;       %lr flip images


%% ---end of user parameters---

addpath('exampledata');
addpath('examplemasks');


date = '20180605';
id = '1603';
flipcmap = 1;

avgfile = 'AVG_20180605_1603.jpg';      %Average projection image
cpfile = 'AVG_20180605_1603_masks.png'; %Cell body segmentation from Cellpose

nstimlevs = 6;
imgwidth = 512;
imgheight = imgwidth;
ImageSize = imgwidth * imgheight;


inpng3 = imread(cpfile);
inpng = inpng3(:,:,1);
avgimg = imread(avgfile);
stdadjust = mat2gray(avgimg);
Low_High = stretchlim(stdadjust,[0.3,0.95]);
stdadjust = imadjust(stdadjust,Low_High,[]);

nsegs = max(max(inpng));      %number of segments
mycells = cell(nsegs,1);      %container to store cell ROIs
ncell = 0;
iscell = zeros(nsegs,1);

for countcells = 1:nsegs
    
    cellROIpix = (inpng == countcells);
    [cellROIy,cellROIx] = find(cellROIpix);
    mycells{countcells}.xpix = cellROIx;
    mycells{countcells}.ypix = cellROIy;
    mycells{countcells}.xcentroid = mean(cellROIx);
    mycells{countcells}.ycentroid = mean(cellROIy);  
    
    if ~isnan(mycells{countcells}.ycentroid)
        iscell(countcells)=1;
        ncell = ncell+1;
    end
    
end

%% cell average response to stimuli
readname = '1603_posq.mat';
load(readname);

minintensity = min(peakvals(:,1:nstimlevs-1),[],2);

nomin = zeros([ImageSize 1]);
denom = zeros([ImageSize 1]);

for i=1:nstimlevs-1
    nomin = nomin + i.* (peakvals(:,i)-minintensity);
    denom = denom + (peakvals(:,i)-minintensity);
end
weightedavg = nomin./denom;
pixelmap = weightedavg;

mycellResponses = cell(1,nsegs);
mycellMaxR = zeros(nsegs,1);
mycellMap = zeros(nsegs,1);
mycellrfs = zeros(nsegs,1);
mycellrfs_accept = zeros(nsegs,1);
mycellfull_accept = zeros(nsegs,1);
mycellr_accept = zeros(nsegs,1);

img_map = zeros(imgwidth,imgheight);
alpha_map = zeros(imgwidth,imgheight);
img_rfs = zeros(imgwidth,imgheight);
alpha_rfs = zeros(imgwidth,imgheight);
alpha_r = zeros(imgwidth,imgheight);

for i=1:nsegs
    ncellpix = length(mycells{i}.xpix);    
    cellxpix = mycells{i}.xpix;
    cellypix = mycells{i}.ypix;
    cellresponse = zeros(1,5);
    cellfull = 0;
    for j=1:ncellpix
        pixindex = sub2ind([imgwidth,imgheight],cellxpix(j),cellypix(j));
        for k=1:nstimlevs-1
            cellresponse(k) = cellresponse(k) + peakvals(pixindex,k);
        end 
        cellfull = cellfull + peakvals(pixindex,nstimlevs); %full field response
    end
    cellresponse = cellresponse./ncellpix;  %normalize
    cellfull = cellfull/ncellpix;
    
    mycellResponses{i} = cellresponse;
    mycellMaxR(i) = max(cellresponse);
    
    if max(cellresponse)>cellrthres
        mycellr_accept(i) = 1;
    end
    
    %rf position   
    maxintensity = max(cellresponse);
    minintensity = min(cellresponse);
    
    cellnomin = 0;
    celldenom = 0;  
    for countstims=1:nstimlevs-1
        cellnomin = cellnomin + countstims.* (cellresponse(countstims)-minintensity);
        celldenom = celldenom + (cellresponse(countstims)-minintensity);
    end
    weightedavg = cellnomin./celldenom;
    mycellMap(i) = weightedavg;
    
    
    %rf sharpness
    rf = weightedavg;
    
    if rf>=2 && rf<=4
        
        %if mycellZ(i)>zthres
            
            dcenter = 0;
            dside = 0;
            ccenter = 0;
            cside = 0;
            for m=1:5
                distance = abs(m-rf);
                if distance<1.5
                    dcenter = dcenter + cellresponse(m);
                    ccenter = ccenter + 1;
                else
                    dside = dside +cellresponse(m);
                    cside = cside + 1;
                end
            end
            
            ratio = (dcenter/ccenter)/(dside/cside);
            mycellrfs(i) = ratio;
            mycellrfs_accept(i) = 1;
        %end

    end
   
    
    %create map images and alpha masks for plotting
    for j=1:ncellpix
        img_map(cellxpix(j),cellypix(j)) = weightedavg;
        alpha_map(cellxpix(j),cellypix(j)) = 1;

        if mycellr_accept(i)
            alpha_r(cellxpix(j),cellypix(j)) = 1;
        end
        
        if mycellrfs_accept(i)
            img_rfs(cellxpix(j),cellypix(j)) = ratio;
            alpha_rfs(cellxpix(j),cellypix(j)) = 1;
        end
        
    end
      
end


threshcells = and(mycellrfs_accept,mycellr_accept);
mycellrfs_clean = mycellrfs(mycellrfs_accept==1);
mycellrfs_thres = mycellrfs(threshcells==1);
ncells_rfs = sum(mycellrfs_accept);
ncells_rfsthres = sum(threshcells);
mrfs = nanmean(mycellrfs_clean);
mrfs_thres = nanmean(mycellrfs_thres);

cellrfs = mycellrfs_clean;
cellrfs_thres = mycellrfs_thres;

%% plotting

%figure 1: show pixel-wise stim position map
img1 = figure(1);

set(img1,'name','Optimal stim position','position',[100 200 512 512],'Color',[1 1 1]);
set(gcf,'units','pixels');
set(gca,'units','pixels');
w_pos = get(gcf, 'position');
set(gca, 'position', [0 0 w_pos(3) w_pos(4)]);

set(img1,'Color',[0 0 0]);
img1.InvertHardcopy = 'off';

posqmap = transpose(reshape(pixelmap,[512,512]));
mapimg = imagesc(posqmap,[1,nstimlevs-1]);
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
if imgxflip
    set(gca, 'XDir','reverse');
end


%figure 2: show cell avg stim position map
img2 = figure(2);
set(img2,'name','Cell average op stim','position',[150 250 512 512],'Color',[0 0 0]);
img2.InvertHardcopy = 'off';
set(gcf,'units','pixels');
w_pos = get(gcf, 'position');
ax1 = axes;
set(ax1,'units','pixels');
set(ax1, 'position', [0 0 w_pos(3) w_pos(4)]);
imagesc(stdadjust);
colormap(ax1,gray);
if imgxflip
    set(gca, 'XDir','reverse');
end

ax2 = axes;
set(ax2,'Color','none','visible','off');
set(ax2,'units','pixels');
set(ax2, 'position', [0 0 w_pos(3) w_pos(4)]);
cellmapimg = imagesc(transpose(img_map),[1 5]);
cellmapimg.AlphaData = transpose(alpha_map);
cellmapimg.AlphaDataMapping = 'scaled';
if flipcmap
    colormap(ax2,flipud(jet));
else
    colormap(ax2,jet);
end
axis off;
if imgxflip
    set(gca, 'XDir','reverse');
end



%figure 3: show cell RF sharpness
img3 = figure(3);
set(img3,'name','Cell RF sharpness','position',[250 150 512 512],'Color',[0 0 0]);
img3.InvertHardcopy = 'off';
set(gcf,'units','pixels');
w_pos = get(gcf, 'position');
ax1 = axes;
set(ax1,'units','pixels');
set(ax1, 'position', [0 0 w_pos(3) w_pos(4)]);
ax1_pos = ax1.Position;
meanimg = imagesc(stdadjust);
colormap(ax1,gray);
if imgxflip
    set(gca, 'XDir','reverse');
end

ax2 = axes;
set(ax2,'Color','none','visible','off');
set(ax2,'units','pixels');
set(ax2, 'position', [0 0 w_pos(3) w_pos(4)]);
cellmapimg = imagesc(transpose(img_rfs));
cellmapimg.AlphaData = transpose(alpha_rfs);
cellmapimg.AlphaDataMapping = 'scaled';
colormap(ax2,hot);
caxis([0.8,1.4]);
axis off;

if imgxflip
    set(gca, 'XDir','reverse');
end



%% display info

disp(['Total segments: ' num2str(nsegs)]);
disp(['Total cells: ' num2str(ncell)]);                         % Total number of cell segments
disp(['Total cells counted: ' num2str(ncells_rfs)]);            % Total number of cells with a central peak stimulus position
disp(['Total cells thresholded: ' num2str(ncells_rfsthres)]);   % Total number of cells with peak response above set threshold
disp(['Mean rf sharpness: ' num2str(mrfs)]);
disp(['Thresholded rf sharpness: ' num2str(mrfs_thres)]);       % Mean rf sharpness after thresholding
