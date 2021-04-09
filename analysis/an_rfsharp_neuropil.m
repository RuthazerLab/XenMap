%{

an_rfsharp_neuropil.m

Calculates receptive field sharpness for neuropil
(RF sharpness = average ?F/F? response to the 2 or 3 stimulus positions
closest to the optimal stimulus position (grid) divided by the average
response to the remaining stimulus positions in the periphery)

%}


%% user parameters

nsamples = 5000;    %number of sample points to show on map
rthres = 2;         %threshold for peak response amplitude (exclude points with peak response below threshold)
nsamples_cdf = 200; %number of sample points to show on cumulative distribution plot

imgxflip = 1;       %lr flip images

%% ---end of user parameters---

addpath('exampledata');
addpath('examplemasks');

nstimlevs = 5;
mapstyle = 'wt';
imgwidth = 512;
imgheight = 512;
ImageSize = imgwidth * imgheight;

maplabel = 'azi';
avgfile = 'AVG_20180605_1603.jpg';
avgimg = imread(avgfile);
stdadjust = mat2gray(avgimg);
%Low_High = stretchlim(stdadjust,[0.3,0.95]);
%stdadjust = imadjust(stdadjust,Low_High,[]);
readname = '1603_posq.mat';
load(readname);

%% find stimulus locations of peak responses
maxintensity = max(max(peakvals));
minintensity = min(peakvals(:,1:nstimlevs-1),[],2);

fullfield = peakvals(:,nstimlevs+1);

if strcmp(mapstyle,'peak')
    [peakvals_next,peakstim] = max(peakvals(:,1:nstimlevs-1),[],2);
    pixelmap = peakstim;
else
    nomin = zeros([ImageSize 1]);
    denom = zeros([ImageSize 1]);
    
    for i=1:nstimlevs-1
        nomin = nomin + i.* (peakvals(:,i)-minintensity);
        denom = denom + (peakvals(:,i)-minintensity);
    end
    weightedavg = nomin./denom;
    pixelmap = weightedavg;
end

%% calculate rf sharpness
means = nan(ImageSize,1);
mcenter = nan(ImageSize,1);
mside = nan(ImageSize,1);

mthres = nan(ImageSize,1);


for k=1:ImageSize
    rf = pixelmap(k);
    meanval = (peakvals(k,1)+peakvals(k,2)+peakvals(k,3)+peakvals(k,4)+peakvals(k,5))/5;
    maxresponse = max(peakvals(k,1:5));
    
    if rf>=2 && rf<=4
        dcenter = 0;
        dside = 0;
        ccenter = 0;
        cside = 0;
        for m=1:5
            distance = abs(m-rf);
            if distance<1.5
                dcenter = dcenter + peakvals(k,m);
                ccenter = ccenter + 1;
            else
                dside = dside + peakvals(k,m);
                cside = cside + 1;
            end
            
        end
        
        mcenter(k) = dcenter/(ccenter*meanval);
        mside(k) = dside/(cside*meanval);
        
        if maxresponse>rthres
            mthres(k) = 1;
        end
        
    end
    
end


%% mask out neuropil & threshold low response values
maskfile = 'AVG_20180605_Exp_1603.tif';
inmask = imread(maskfile);
fulltp = transpose(reshape(fullfield,[512,512]));
fullmasked = fulltp;
fullmasked(inmask==0) = NaN;
meanfull = nanmean(reshape(fullmasked,[1,ImageSize]));
stdfull = nanstd(reshape(fullmasked,[1,ImageSize]));

mcentermasked = transpose(reshape(mcenter,[512,512]));
msidemasked = transpose(reshape(mside,[512,512]));
mcentermasked(inmask==0)=NaN;
msidemasked(inmask==0)=NaN;

minds = find(~isnan(mcentermasked));
minds_thres = find(and(~isnan(mcentermasked),~isnan(reshape(mthres,[512,512]))));

mcenterflat = reshape(mcentermasked,[1,ImageSize]);
msideflat = reshape(msidemasked,[1,ImageSize]);

mcenterflatthres = mcenterflat(~isnan(mthres));
msideflatthres = msideflat(~isnan(mthres));

mcenterflat = mcenterflat(~isnan(mcenterflat));
msideflat = msideflat(~isnan(msideflat));
mcenterflatthres = mcenterflatthres(~isnan(mcenterflatthres));
msideflatthres = msideflatthres(~isnan(msideflatthres));

onesvector = ones(1,length(mcenterflat));
twosvector = 2*ones(1,length(msideflat));

mean1 = mean(mcenterflat);
mean2 = mean(msideflat);
meanplot = [mean1,mean2];
stderr1 = std(mcenterflat)/sqrt(length(mcenterflat));
stderr2 = std(msideflat)/sqrt(length(msideflat));
errorplot = [stderr1,stderr2];

mratio = mcenterflat./msideflat;
rfsmean = mean(mratio);
rfsmedian = median(mratio);

mratio_thres = mcenterflatthres./msideflatthres;
rfsmeanthres = mean(mratio_thres);
rfsmedianthres = median(mratio_thres);

nraw = length(mratio);
nthres = length(mratio_thres);

rfssample = datasample(mratio,nsamples_cdf,'Replace',false);
rfsthressample = datasample(mratio_thres,nsamples_cdf,'Replace',false);


%% display plot

npoints = min(length(minds_thres),nsamples);
samples = randsample(length(minds_thres),npoints);
sampleinds = minds_thres(samples);
sample_mratio = mratio_thres(samples);

img1 = figure(1);
set(img1,'name','RF sharpness','position',[250 150 512 512],'Color',[0 0 0]);
img1.InvertHardcopy = 'off';
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
[sampleindsy,sampleindsx] = ind2sub([512,512],sampleinds);
scatter(sampleindsx,sampleindsy,8,sample_mratio,'filled');
axis off;
set(gca, 'YDir','reverse');
if imgxflip
    set(gca, 'XDir','reverse');
end
xlim([0,512]);
ylim([0,512]);
caxis([0.8,1.4]);
colormap(ax2,hot);

img2 = figure(2);
set(img2,'name','RF sharpness: cumulative distribution')
cdfplot(rfsthressample);
xlim([0.6,2.2]);
title('RF sharpness: cumulative distribution');
