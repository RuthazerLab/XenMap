%{

an_discont.m

Calculates local discontinuity: the mean difference in receptive field position 
(phase) of a pixel to all neighboring pixels within a defined pixel radius

1. Mask input phase map with neuropil mask and threshold for SNR
2. Randomly select pixels from neuropil area to evaluate
3. Generate neighborhod circle (poly2mask->bwtrace)
4. Acquire coordinates of points within circle
5. Calculate and average phase distance
6. Display scatter plot of points colored by discontinuity values

%}


%% user parameters

radpixs = 7;            %neighborhood radius (in pixels)
ncirclepoints = 20;     %number of points on sampling circle circumference
snrthres = 1 ;          %SNR threshold
nneighbors = 20;        %least number of neighboring pixels (exclude pixels with insufficient #neighbors after thresholding for SNR)
nsamples = 5000;      	%number of pixels to sample

imgxflip = 1;       %lr flip images

%% ---end of user parameters---

addpath('exampledata');
addpath('examplemasks');

imgfile = 'AVG_Exp_1732_slice2.jpg';    %average projection image
maskfile = 'AVG_Exp_1732_slice2.tif';   %neuropil mask
mapfile = '20180722_1727-1732_sl2_filt=1.mat';  %phase map
datafile = 'StimulusData_20180722_17.27.mat'; %experiment info

load(mapfile);
load(datafile);
pblank = setup.DriftSettings.Blanktime/(setup.DriftSettings.Blanktime+1/setup.DriftSettings.Frequency);
orientation = setup.DriftSettings.Orientation;
if strcmp(orientation,'horizontal')
    clims = [-pi -pi+2*pi*(1-pblank)*0.5];
    axislim = 85;
    stimaxis = 'elevation';
    ylims = [-pi,-pi+2*pi*(1-pblank)*0.5];
else
    clims = [-pi -pi+2*pi*(1-pblank)];
    axislim = 110;
    stimaxis = 'azimuth';
    ylims = [-pi,-pi+2*pi*(1-pblank)];
end

%% apply neuropil mask and threshold by SNR

snrmask_thres = zeros(512);
snrmask_thres(snrmap>=snrthres)=1;
snrmask_thres = transpose(snrmask_thres);
maskedmap = map_sub;
maskedmap(snrmask_thres==0)=0;

%read mask from file
inmask = imread(maskfile);
maskedmap(inmask==0)=0;

%calculate average neuropil SNR
snr_masked = snrmap;
snr_masked(inmask==0) = NaN;
snr_masked(snrmap<snrthres)=NaN; %%2021.1.15: calculate mean SNR after threshold
snr_mean = nanmean(nanmean(snr_masked));
disp(['Average neuropil SNR is ' num2str(snr_mean)]);

%randomly sample n points and create histograms of tectal distance vs vf distance
maskpoints = find(maskedmap~=0);
npoints = min(length(maskpoints),nsamples);
sampleinds = randsample(maskpoints,npoints);
t = linspace(0, 2*pi, ncirclepoints);

vfdistsmain = {length(radpixs)};

r = radpixs;
%vfdists = NaN(nsamples,round(2*pi*(radpixs(crtrad)+1)));
vfdists = NaN(1,npoints);
nkeptpoints = 0;
for crtpoint = 1:npoints
    crtpointphase = maskedmap(sampleinds(crtpoint));
    [crty,crtx]=ind2sub([512,512],sampleinds(crtpoint));
    c = [crtx,crty];
    % get circular mask
    BW = poly2mask(r*cos(t)+c(1), r*sin(t)+c(2), 512,512);
    neighpoints = find(BW~=0);
    vfdist_avg = 0;
    vfdist_count = 0;
    for neighc = 1:length(neighpoints)
        comparephase = maskedmap(neighpoints(neighc));
        if comparephase~=0
            phasediff = abs(comparephase-crtpointphase);
            vfdist = axislim*phasediff/(2*pi*(1-pblank));
            vfdist_avg = vfdist_avg + vfdist;
            vfdist_count = vfdist_count+1;
        end
    end
    if vfdist_count>nneighbors
        vfdist_avg = vfdist_avg/vfdist_count;
        vfdists(crtpoint) = vfdist_avg;
        nkeptpoints = nkeptpoints+1;
    end
    
end
vfdists_avg = nanmean(vfdists);
vfdists_med = nanmedian(vfdists);
%disp(['Average neuropil RF distance is ' num2str(vfdists_avg)]);

[sampleindsy,sampleindsx] = ind2sub([512,512],sampleinds);

img1 = figure(1);
set(img1,'position',[100 200 700 500]);

subplot(2,2,[3,4]);
histogram(vfdists);
title({['Mean Discontinuity = ' num2str(vfdists_avg) '°  Mean SNR=' num2str(snr_mean) ]});
xlim([0,0.2*axislim]);
ylabel('Discontinuity');

subtightplot(2,2,1,0.1);
imagesc(maskedmap);
if imgxflip
    set(gca, 'XDir','reverse');
end
axis off;
colorbar;
caxis(clims);
title('Phase map')


subtightplot(2,2,2,0.1)
scatter(sampleindsx,sampleindsy,2,vfdists,'filled');
set(gca, 'YDir','reverse');
if imgxflip
    set(gca, 'XDir','reverse');
end
title('Discontinuity distribution');
xlim([0,512]);
ylim([0,512]);
colorbar;
caxis([0,0.15*axislim]);
colormap jet



disp (['Mean neuropil discontinuity is ' num2str(vfdists_avg) '°']);



