close all; clear; clc;

%%%%%%%%%%360LI algorithm, benchtop, PCD operated as if it only had a
%%%%%%%%%%single row

%Initialize parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
detz = 4.3; %detector z-coverage
pitch = 1; % for 360 degree rotation
omega = 4; %deg/s
totframes = 4500;
frames = 900; %frames per 360 degrees
numSweeps = totframes/frames; %scan range = numSweeps*v
slices = 10; %recon slices per detz
PixelBinning = 4;

v = detz*pitch; %obj z-translation (mm) per 360 degree rotation
%%
%Import data and generate initial sinogram with corrections %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
objFilename = 'data\20220331_Helical\1_TE.EVI';
%objFilename = 'data\20220330_Helical\1_TE.raw';
airFilename = 'data\20220331_Helical\air_TE.EVI';
%airFilename = 'data\20220330_Helical\air_TE.raw';

obj = MgReadEviData(objFilename);
%obj = MgReadRawFile(objFilename,64,5120,totframes,0,0,'float32');
air = MgReadEviData(airFilename);
%air = MgReadRawFile(airFilename,64,5120,910,0,0,'uint32');

air = mean(air, 3);

temp = air./obj;
clear air
clear obj

sgm = log(temp);
clear temp

sgm(isnan(sgm)) = 0;
sgm(isinf(sgm)) = 0;
%%
%Panel correction
load('coefficients.mat');
sgm = MgPolyval(cali_PMMA_TE, sgm);
sgm(isnan(sgm)) = 0;
sgm(isinf(sgm)) = 0;
%%%%%

%Panel gap interpolation
sgm = MgInterpSinogram(sgm, 2, 2, 4);

sgm = MgResliceProjectionToSinogram(sgm, 3, 60, 1);%I module PCD 1 slice
sgm = MgRebinSinogram(sgm,PixelBinning); %pixel binning

%% 4/1 360LI redo from scratch

sgmHelical=zeros(frames,5120/PixelBinning,slices*pitch*numSweeps);
zProjs=0:v/frames:numSweeps*v-v/frames;
thetaProjs = 0:(360/frames):numSweeps*360-360/frames;

for i=1:(slices*pitch*numSweeps-1)
    z=(i-1)*detz/slices;
    if z<v % extrap %%%%%%%%%%%%%%%%change logical range
    elseif z>=(numSweeps-1)*v %extrap %%%%%%%%%%%%%%%% change logical range
    else %interpolation z-range
        for j=1:frames
            sweepNumber = floor(z/v);
            if z<(sweepNumber*v+j*v/frames)
                sweepNumber=sweepNumber-1;
            end
            projIndexLow = j+(sweepNumber)*frames;
            projIndexHigh = projIndexLow+frames;
            zLow = zProjs(projIndexLow);
            zHigh = zProjs(projIndexHigh);
            weightLow = (zHigh-z)/v;
            weightHigh = (z-zLow)/v;
            sgmHelical(j,:,i) = weightLow.*sgm(projIndexLow,:)+weightHigh.*sgm(projIndexHigh,:);
        end
    end
end

%%
filename = sprintf('sgm/20220331_Helical/sgm_1_%d-%d-%d.raw', size(sgmHelical, 2), size(sgmHelical, 1), size(sgmHelical,3));
MgSaveRawFile(filename, sgmHelical);

%%
%!mgfbp config_mgfbp_Treb_I.jsonc

% % close all; clear; clc;
% img = MgReadRawFile('img/20220324_QA_center_narrow/img_bh_1_512-512-60.raw', 512, 512, 60);
% for i=1:60
%     imgtemp = img(:,:,i);
%     imgtemp = MgRingCorrection(imgtemp, -1000000, 1000000, 100, 20, 70);
%     img(:,:,i) = imgtemp;
% end
% 
% MgSaveRawFile('img/20220324_QA_center_narrow/img_bh_ring_1_512-512-60.raw', img);
