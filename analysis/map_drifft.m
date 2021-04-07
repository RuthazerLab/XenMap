%{

map_drifft.m

For calculating phase maps: Performs fast fourier transform on
pixel-wise calcium traces and calculates phase at stimulus frequency

Input: raw movie (16 bit unsigned, little endian byte order)
(response to repeated presentations of drifting bar stimulus)

Output: .mat file with 4 data fields
-phasemap: map of phase at stimulus freq
-lummask: average projection image of first 500 frames
-snrmap: signal to noise ratio map
-magmap: map of response amplitude at stimulus freq

%}

%% user parameters

gaussfilt = 1;      %sigma value for Gaussian filter (if input 0, then skip filtering step)

filename = '';      %path to raw movie
imgwidth = 512;     %size of raw image
imgheight = 512;
totstimframes = 1000;   %total number of frames to read from movie
CaptureRate = 6;    %frame rate of movie
StimRate = 0.05;     %stimulus presentation frequency

diffbin = 5;        %bin size for calculating first differential
dist = 3;           %frame distance between bins subtracted for first differential

writename = '';     %name of output data file

%% ---end of user parameters---

fid = fopen(filename,'r','l');
wholeclip = fread(fid,[1 ImageSize*totstimframes],'uint16','ieee-le');
wholeclip = reshape(wholeclip,[ImageSize,totstimframes]);

%Gauss filter
if gaussfilt
    for i = 1:totstimframes
        filtframe = reshape(wholeclip(:,i),[imgwidth,imgheight]);
        filtframe= imgaussfilt(filtframe,gaussfilt);
        wholeclip(:,i) = reshape(filtframe,[1,ImageSize]);
    end
end


%Generate average projection of first 500 frames
lummask = wholeclip(ImageSize,500);
lummask = mean(lummask,2);
lummask = transpose(reshape(lummask,[imgwidth,imgheight]));


%Calculate first differential traces

num_imgsn = totstimframes-(diffbin-1);
all_imgsn = zeros(ImageSize,num_imgsn);

for ii = 1:diffbin
    all_imgsn = all_imgsn + wholeclip(:,ii:num_imgsn-1+ii);
end

all_imgsn = all_imgsn./diffbin;
difff = all_imgsn(:,dist+2:num_imgsn) - all_imgsn(:,1:num_imgsn-dist-1);
med_difff = median(difff);
multiplier = ones(ImageSize,1);
difff = difff - multiplier*med_difff;
totstimframes = num_imgsn-dist-1;


%Calculate phase map
f = ((0:1/totstimframes:1-1/totstimframes)*CaptureRate).';
[~,freqind] = min(abs(f-StimRate));
[~,harm2] = min(abs(f-StimRate*2));
[~,harm3] = min(abs(f-StimRate*3));
[~,harm4] = min(abs(f-StimRate*4));
disp('Calculating phase map...');
phasemap = zeros(ImageSize,1);
snrmap = zeros(ImageSize,1);
magmap = zeros(ImageSize,1);

for i=1:ImageSize
    
    y = fft(difff(i,:));
    
    
    phasemap(i) = -1*unwrap(angle(y(freqind)));
    %Calculate SNR: A/SD where A is amplitude at stim freq, SD is stdev of
    %amplitudes above stim freq (excluding first 3 harmonics)
    magnitudes = abs(y);
    stimmag = magnitudes(freqind);
    magnitudes(harm4)=[];
    magnitudes(harm3)=[];
    magnitudes(harm2)=[];
    magnitudes(1:freqind)=[];
    magmap(i) = stimmag;
    snrmap(i) = stimmag/std(magnitudes);
end

save(writename,'phasemap','lummask','snrmap','magmap');

fclose(fid);