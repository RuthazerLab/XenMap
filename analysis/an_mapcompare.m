%{

an_mapcompare.m

Calculates difference score between 2 maps (or a map and its scrambled
version)

%}


%% user parameters

gaussfilt = 3;         %perform Gaussian filter on both maps before comparing
dobootstrap = 1000;    %perform bootstrap with scrambled image (set #iterations)

imgxflip = 1;       %lr flip images

%% ---end of user parameters---

addpath('exampledata');
addpath('examplemasks');
addpath('utilities');

stimaxis = 'azimuth';

if strcmp(stimaxis,'elevation')
    clims = [-pi -pi+2*pi*0.5*0.5];
else
    clims = [-pi -pi+2*pi*0.5];
end

targetmapfile = '20180722_1727-1732_sl2_filt=1.mat';
load(targetmapfile);
tgtmap = map_sub;
if gaussfilt
    tgtmap = imgaussfilt(tgtmap,gaussfilt);
end
tgtsnr = transpose(reshape(snrmap,[512,512]));
tgtlum = lummask2;

compmapfile = '1811_posq.mat';
compmaskfile = '1811_posqmask.mat';
load(compmapfile);
load(compmaskfile);
cmpmap = transpose(reshape(pixelmap,[512,512]));
if gaussfilt
    cmpmap = imgaussfilt(cmpmap,gaussfilt);
end
cmplum = fullmask;

%display map 1: phase map
img1 = drawMap(1,map_sub,snrmap,clims,[512 512]);
set(img1,'name','Phase map','position',[100 200 512 512]);
if imgxflip
    set(gca, 'XDir','reverse');
end

%display map 2: grid map
img2 = drawMapFlash(2,pixelmap,cmplum,flipcmap,clims,[512 512]);
set(img2,'name','Grid map','position',[150 250 512 512]);
if imgxflip
    set(gca, 'XDir','reverse');
end

%adjust grid map values to match phase map range
cmpmap = 6-cmpmap;
cmpmap = rescale(cmpmap,min(clims),max(clims));

%segment neuropil from mask
maskfile = 'AVG_Exp_1732_slice2.tif';
inmask = imread(maskfile);
bwmask = inmask>0;
tgtmap(bwmask==0) = NaN;

% register and transform maps for comparison
[optimizer,metric] = imregconfig('monomodal');
tform1 = imregtform(tgtlum, cmplum,'translation',optimizer,metric);
Rfixed = imref2d(size(cmplum));
cmpmask = imwarp(bwmask,tform1,'OutputView',Rfixed);

lummaskshow = mat2gray(tgtlum);
Low_High = stretchlim(lummaskshow,[0.3,0.9]);
lummaskshow = imadjust(lummaskshow,Low_High,[]);

cmpmapshow = mat2gray(cmplum);
Low_High = stretchlim(cmpmapshow,[0.3,0.9]);
cmpmapshow = imadjust(cmpmapshow,Low_High,[]);

% img3 = figure(3);
% set(img3,'name','Registration','position',[200 300 800 300],'Color',[1 1 1]);
% subplot(1,2,1)
% regimg = imshowpair(cmpmapshow, lummaskshow, 'Scaling', 'joint');
% title('Registration');
% subplot(1,2,2)
% regimg2 = imshowpair(cmpmapshow, lummaskshow, 'montage');

%compute difference score for phase vs grid map
cmpmap(cmpmask==0) = NaN;
tgtmap(isnan(cmpmap))=NaN;
cmpmap(isnan(tgtmap))=NaN;
diff1 = abs(tgtmap - cmpmap);
diff1score = nanmean(nanmean(diff1));

%produce shuffled version of phase map
targetpoints = find(~isnan(cmpmap));
shuffledpoints = targetpoints(randperm(length(targetpoints)));
shuffledmap = nan([512,512]);
shuffledvals = cmpmap(shuffledpoints);
shuffledmap(targetpoints)=shuffledvals;

diff3 = abs(tgtmap - shuffledmap);
diff3score = nanmean(nanmean(diff3));


 %bootstrap
 if dobootstrap
     bootscores = nan(1,dobootstrap);
     for i=1:dobootstrap
         shuffledpoints2 = targetpoints(randperm(length(targetpoints)));
         shuffledmap2 = nan([512,512]);
         shuffledvals2 = tgtmap(shuffledpoints2);
         shuffledmap2(targetpoints)=shuffledvals2;
         diff4 = abs(tgtmap - shuffledmap2);
         bootscores(i) = nanmean(nanmean(diff4));
     end
 end
 
 %display comparison results
 if dobootstrap
     img4 = figure(4);
     histogram(bootscores);
     title(['Scrambled difference scores - ' num2str(dobootstrap) ' bootstrap iterations']);
     set(img4,'position',[400 100 800 300]);
 end

 
 
 img5 = figure(5);
 set(img5,'position',[300 200 500 330]);
 
 subtightplot(2,3,1);
 imagesc(tgtmap);
 title('Target');
 colormap jet;
 caxis(clims);
if imgxflip
    set(gca, 'XDir','reverse');
end
 text(500,480,'Phase map','Color','w');
 axis off;
 
 subtightplot(2,3,2);
 imagesc(cmpmap);
 title('Compare');
 caxis(clims);
if imgxflip
    set(gca, 'XDir','reverse');
end
 text(500,480,'Grid map','Color','w');
 axis off;
 
 subtightplot(2,3,3);
 imagesc(diff1);
 title('Difference');
 caxis([-0.2 1.5]);
if imgxflip
    set(gca, 'XDir','reverse');
end
 text(500,480,['Diff score: ' num2str(diff1score)],'Color','w');
 axis off;
 
 subtightplot(2,3,4);
 imagesc(tgtmap);
 colormap jet;
 caxis(clims);
if imgxflip
    set(gca, 'XDir','reverse');
end
 text(500,480,'Phase map','Color','w');
 axis off;

 subtightplot(2,3,5);
 imagesc(shuffledmap);
 caxis(clims);
if imgxflip
    set(gca, 'XDir','reverse');
end
 text(500,480,'Scrambled Phase','Color','w');
 axis off;
 
 subtightplot(2,3,6);
 imagesc(diff3);
 caxis([-0.2 1.5]);
if imgxflip
    set(gca, 'XDir','reverse');
end
 text(500,480,['Diff score: ' num2str(diff3score)],'Color','w');
 axis off;



 
 
 