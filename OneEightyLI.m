close all; clear; clc;

%%%%%%%%%%180LI algorithm, benchtop, PCD operated as if it only had a
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

%extra 180 LI parameters
DetectorOffcenter = 43;
SDD = 1041.4;

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

%% 4/4 180LI det and gamma negative and no truncation
% notes, plus or minus DetectorOffcenter?
% notes, gamma plus or minus?
% notes, extrap z conditions for if statements

sgmHelical=zeros(frames,5120/PixelBinning,slices*pitch*numSweeps);
zProjs=0:v/frames:numSweeps*v-v/frames;
thetaProjs = 0:(360/frames):numSweeps*360-360/frames;

for i=1:(slices*pitch*numSweeps-1) % through recon z
    z=(i-1)*detz/slices;
    if z<v % extrap
    elseif z>((numSweeps-1)*v) %extrap
    else %interpolation z-range
        sweepNumber = floor(z/v);
        for j=1:frames % through beta
            if z<(sweepNumber*v+j*v/frames)
                sweepNumber=sweepNumber-1;
            end
            projIndexLow = j+(sweepNumber)*frames; % initial beta
            zLow = zProjs(projIndexLow);

            for k=(2*DetectorOffcenter/(0.1*PixelBinning)+1):(5120/PixelBinning-2*DetectorOffcenter/(0.1*PixelBinning)-1) %through pixel index with conjugate rays (transformed to gamma)
                gamma = atand(((k-0.5)*0.1*PixelBinning-5120*0.1/2+DetectorOffcenter)/SDD); % gamma in degrees
            
                projIndexHigh = round(projIndexLow+frames/2+2*gamma*frames/360); % beta plus 180 degrees plus 2*gamma (rounded)
                zHigh = zProjs(projIndexHigh);
                
                if z>zHigh
                    projIndexHighnew = projIndexLow+frames;
                    projIndexLownew = projIndexHigh;
                    zHighnew = zLow+v;
                    zLownew = zHigh;
                    weightLownew = (zHighnew-z)/(zHighnew-zLownew);
                    weightHighnew = (z-zLownew)/(zHighnew-zLownew);
                    kLownew = (-((k-0.5)*0.1*PixelBinning-5120*0.1/2+DetectorOffcenter)-DetectorOffcenter+5120*0.1/2)/(0.1*PixelBinning)+0.5;
                    sgmHelical(j,k,i) = weightLownew.*sgm(projIndexLownew,round(kLownew))+weightHighnew.*sgm(projIndexHighnew,round(k));
                else      
                    weightLow = (zHigh-z)/(zHigh-zLow);
                    weightHigh = (z-zLow)/(zHigh-zLow);

                    %kHigh=tand(-gamma)*SDD/(0.1*PixelBinning);
                    kHigh = (-((k-0.5)*0.1*PixelBinning-5120*0.1/2+DetectorOffcenter)-DetectorOffcenter+5120*0.1/2)/(0.1*PixelBinning)+0.5;
                    sgmHelical(j,k,i) = weightLow.*sgm(projIndexLow,k)+weightHigh.*sgm(projIndexHigh,round(kHigh));
                end
            end
            for k=1:(2*DetectorOffcenter/(0.1*PixelBinning)) % 360 LI when no conjugate rays left side det
                projIndexHigh = projIndexLow+frames;
                zHigh = zProjs(projIndexHigh);
                weightLow = (zHigh-z)/v;
                weightHigh = (z-zLow)/v;
                sgmHelical(j,k,i) = weightLow.*sgm(projIndexLow,k)+weightHigh.*sgm(projIndexHigh,k);
            end
            for k=(5120/PixelBinning-2*DetectorOffcenter/(0.1*PixelBinning)):5120/PixelBinning % 360 LI when no conjugate rays right side det
                projIndexHigh = projIndexLow+frames;
                zHigh = zProjs(projIndexHigh);
                weightLow = (zHigh-z)/v;
                weightHigh = (z-zLow)/v;
                sgmHelical(j,k,i) = weightLow.*sgm(projIndexLow,k)+weightHigh.*sgm(projIndexHigh,k);
            end
        end
    end
end

%%
filename = sprintf('sgm/20220331_Helical/Det_Gamma_neg_180LIsgm_1_%d-%d-%d.raw', size(sgmHelical, 2), size(sgmHelical, 1), size(sgmHelical,3));
MgSaveRawFile(filename, sgmHelical);

%%
!mgfbp config_mgfbp_Treb_I_ReverseHelical.jsonc

% close all; clear; clc;
img = MgReadRawFile('img/ReverseHelical_tube_on/360LIimg_bh_1_512-512-48.raw', 512, 512, 48);
for i=1:48
    imgtemp = img(:,:,i);
    imgtemp = MgRingCorrection(imgtemp, -0.8, 100, 0.05, 20, 70);
    img(:,:,i) = imgtemp;
end

MgSaveRawFile('img/ReverseHelical_tube_on/360LIimg_bh_ring_1_512-512-48.raw', img);

