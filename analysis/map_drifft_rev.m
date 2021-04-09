%{

map_drifft_rev.m

For calculating phase maps: Calculates absolute phase map by taking 
difference of phase maps acquired from opposite direction drifting bar stimuli

Input: pair of .mat files containing relative phase maps to be corrected 
(extracted from opposite direction drifting bar experiments)


*For registering two experiments, option to use either Matlab library
method for rigid registration or NoRMCorre package for non-rigid registration
github.com/flatironinstitute/NoRMCorre

%}


%% user parameters

imgwidth=512;   %input image size
imgheight=512;

stimaxis = 'azimuth'; %stimulus axis (azimuth or elevation)
driftTime = 10;     %time drifting bar stimulus is on screen during each cycle
driftBlank = 10;    %time drifting bar stimulus is off screen during each cycle

register = 1;   %0=none / 1=use Matlab default / 2=use nonrigid registration from NoRMCorre / 3=use rigid registration from NoRMCorre

imgxflip = 1;       %lr flip images

%% ---end of user parameters---

addpath('exampledata');
addpath('utilities');

showreg = 0;

map1id = '1727_diff_phasemap_sl2_filt=1.mat';
map2id = '1732_diff_phasemap_sl2_filt=1.mat';

load(map1id);
phasemap1 = phasemap;
lummask1 = lummask;
load(map2id);
phasemap2 = phasemap;
lummask2 = lummask;

pblank = driftBlank/(driftTime + driftBlank);

if strcmp(stimaxis, 'elevation')
    clims = [-pi -pi+2*pi*(1-pblank)*0.5];
else
    clims = [-pi -pi+2*pi*(1-pblank)];
end


phasemap1 = transpose(reshape(phasemap1,[imgwidth,imgheight]));
phasemap2 = transpose(reshape(phasemap2,[imgwidth,imgheight]));


%% register 2 maps based on average projection images
if register ==1
    [optimizer,metric] = imregconfig('monomodal');
    tformSimilarity = imregtform(lummask1, lummask2,'translation',optimizer,metric);
    Rfixed = imref2d(size(lummask2));
    movingRegisteredRigid = imwarp(lummask1,tformSimilarity,'OutputView',Rfixed);
    
    lummaskshow = mat2gray(lummask2);
    Low_High = stretchlim(lummaskshow,[0.3,0.9]);
    lummaskshow = imadjust(lummaskshow,Low_High,[]);
    
    movingRegisteredRigid = mat2gray(movingRegisteredRigid);
    Low_High = stretchlim(movingRegisteredRigid,[0.3,0.9]);
    movingRegisteredRigid = imadjust(movingRegisteredRigid,Low_High,[]);
    
    if showreg
        regwin = figure(1);
        regimg = imshowpair(movingRegisteredRigid, lummaskshow, 'Scaling', 'joint');
        title('Registration');
    end
    
    mapimg1 = movingRegisteredRigid;
    mapimg2 = lummaskshow;
    
    phasemap1 = imwarp(phasemap1,tformSimilarity,'OutputView',Rfixed);
    
    %code using NoRMCorr
elseif register==2 || register==3
    if register == 2
        options = NoRMCorreSetParms('d1',size(lummask1,1),'d2',size(lummask1,2),'grid_size',[64,64],'overlap_pre',[16,16],'mot_uf',4,'bin_width',50,'max_shift',[40,40],'max_dev',[8,8],'us_fac',50,'init_batch',200,'shifts_method','FFT','correct_bidir',false);
        [lummask1_reg,shifts1,template,options_update] = normcorre(single(lummask1),options,single(lummask2));
    else
        options = NoRMCorreSetParms('d1',size(lummask1,1),'d2',size(lummask1,2),'bin_width',50,'max_shift',40,'us_fac',50,'init_batch',200,'correct_bidir',false,'shifts_method','FFT');
        [lummask1_reg,shifts1,template,options_update] = normcorre(single(lummask1),options,single(lummask2));
    end
    
    lummaskshow = mat2gray(lummask2);
    Low_High = stretchlim(lummaskshow,[0.3,0.9]);
    lummaskshow = imadjust(lummaskshow,Low_High,[]);
    
    lummask1_reg = mat2gray(lummask1_reg);
    Low_High = stretchlim(lummask1_reg,[0.3,0.9]);
    lummask1_reg = imadjust(lummask1_reg,Low_High,[]);
    
    if showreg
        regwin = figure(1);
        set(regwin,'name','Registration','position',[500 200 800 300],'Color',[1 1 1]);
        subplot(1,2,1)
        regimg = imshowpair(lummask1_reg, lummaskshow, 'Scaling', 'joint');
        title('Registration');
        subplot(1,2,2)
        regimg2 = imshowpair(lummask1_reg, lummaskshow, 'montage');
    end
    
    mapimg1 = mat2gray(lummask1);
    Low_High = stretchlim(mapimg1,[0.3,0.9]);
    mapimg1 = imadjust(mapimg1,Low_High,[]);
    
    mapimg2 = lummaskshow;
    
    phasemap1 = apply_shifts(single(phasemap1),shifts1,options_update);
    
end


%% calculate absolute phse map and display
map_sub = (phasemap1 - phasemap2 - pblank*2*pi) ./2;


img2 = drawMap(2,map_sub,snrmap,clims,[512 512]);
set(img2,'name','Phase map','position',[100 200 512 512]);
if imgxflip
    set(gca, 'XDir','reverse');
end