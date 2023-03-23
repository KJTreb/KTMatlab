close all; clear; clc;

%%%%%%%%%%360LI algorithm for reverse extended helical trajectory, benchtop, PCD operated as if it only had a
%%%%%%%%%%single row

%Initialize parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
detz = 4.3; %detector z-coverage
%pitch = 0.6; % for shortscan plus rest
fps = 10;
omega = 4.16; %deg/s
totFrames = 4350;
scanFrames = 500; %frames per short scan
restFrames = 50; %frames per rest period
numSweeps = 8; %scan range = numSweeps*v
slices = 10; %recon slices per detz
PixelBinning = 4;
zPerFrame = 0.001875*25.4/fps; % mm/frame

zScan = scanFrames*zPerFrame; %obj z-translation (mm) per short scan
zRest = restFrames*zPerFrame; %obj z-translation (mm) per rest period

%extra 180 LI parameters
DetectorOffcenter = 43;
SDD = 1041.4;

%%
%Import data and generate initial sinogram with corrections %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
objFilename = 'data\20220405_ReverseHelical\1_TE.EVI';
airFilename = 'data\20220331_Helical\air_TE.EVI';

obj = MgReadEviData(objFilename);
air = MgReadEviData(airFilename);

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

sgmOG=sgm;
%% 4/18 360LI reverse helical trajectory
% notes, extrap z conditions for if statements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%known delay between data acquisition and start of phantom rotation
delay=50; %(50) can be up to 50, if higher then reaquire
sgmStageTwo=sgmOG((delay+1):(4400-(50-delay)),:);

% empirical offset between start of rotation and start of fram acquisition
% (up to 0.5 times frame time)
sgm = zeros(size(sgmStageTwo));
offset = 0.1; % XXX/fps (empirical, between 0 and 1, used as weighting factor)
for p=1:(length(sgm(:,1))-1)
    sgm(p,:) = (1-offset)*sgmStageTwo(p,:) + offset*sgmStageTwo(p+1,:);
end
%note don't need to account for reverse sweep direction because we are
%working with the original full-acquisition sinogram prior to "flipping"
%the direction in the construction of the final sinogram used for recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sgmHelical=zeros(scanFrames,5120/PixelBinning,floor(zPerFrame*totFrames*slices/detz));
zProjs=0:zPerFrame:zPerFrame*(totFrames-1);

thetaProjs1 = 0:(208/scanFrames):(208/scanFrames*(scanFrames-1));
thetaProjs2 = [thetaProjs1,208/scanFrames*(scanFrames-1)*ones(1,50)];
thetaProjs3 = (208/scanFrames*(scanFrames-1)):(-208/scanFrames):0;
thetaProjs4 = cat(2, thetaProjs2, thetaProjs3, zeros(1,50));

thetaProjs = cat(2, thetaProjs4, thetaProjs4, thetaProjs4, thetaProjs2, thetaProjs3); %in degrees, source angle between 0 and ~208

for i=7:37%1:(floor(zPerFrame*totFrames*slices/detz)-1) % through recon z
    z=(i-1)*detz/slices;
    if z<(zScan+zRest) % extrap
    elseif z>((numSweeps-1)*(zScan+zRest)) %extrap
    else %interpolation z-range
        sweepNumber = floor(z/(zScan+zRest));
        for j=1:scanFrames % through beta 
            
            if rem(sweepNumber,2)==1 %if reverse sweep, then reverse direction of j
                frameIndex=scanFrames-j+1;
            else
                frameIndex=j;
            end
            
            if z<(sweepNumber*(zScan+zRest)+frameIndex*zPerFrame) % check if recon z is still below chosen frame's z from original sinogram
                sweepNumber=sweepNumber-1;
                frameIndex=scanFrames-frameIndex+1;
            end
            
            if j==1
                if (4.9959<z && z<5.2387) %inside rest zone and at angle 0
                    projIndexLow = find(zProjs>4.9959 & zProjs<5.2387); %Bass Ackwards way
                    for k=1:length(projIndexLow)
                        sgmHelical(j,:,i) = sgmHelical(j,:,i)+sgm(projIndexLow(k),:);
                    end
                elseif (10.2346<z && z<10.4775)
                    projIndexLow = find(zProjs>10.2346 & zProjs<10.4775); %Bass Ackwards way
                    for k=1:length(projIndexLow)
                        sgmHelical(j,:,i) = sgmHelical(j,:,i)+sgm(projIndexLow(k),:);
                    end
                elseif (15.4734<z && z<15.7162)
                    projIndexLow = find(zProjs>15.4734 & zProjs<15.7162); %Bass Ackwards way
                    for k=1:length(projIndexLow)
                        sgmHelical(j,:,i) = sgmHelical(j,:,i)+sgm(projIndexLow(k),:);
                    end
                else
                    projIndexLow = frameIndex+(sweepNumber)*(scanFrames+restFrames); % initial beta
                    zLow = zProjs(projIndexLow);
                    %360LI
                    projIndexHigh = projIndexLow+restFrames+2*(scanFrames-frameIndex)+1; % minus rest frames ???
                    zHigh = zProjs(projIndexHigh);
                    
                    if zHigh < z % brute force and ignorance method of making sure its interpolation and not extrapolation
                        sweepNumber=sweepNumber+1;
                        frameIndex=scanFrames-frameIndex+1;
                        projIndexLow = frameIndex+(sweepNumber)*(scanFrames+restFrames); % initial beta
                        zLow = zProjs(projIndexLow);
                        projIndexHigh = projIndexLow+restFrames+2*(scanFrames-frameIndex)+1;
                        zHigh = zProjs(projIndexHigh);
                    end
%%%
                    if (zLow>4.9959 && zLow<10.2346)
                        projIndexLow = find(zProjs>4.9959 & zProjs<5.2387); %Bass Ackwards way
                        projIndexHigh = find(zProjs>10.2346 & zProjs<10.4775);
                        weightLow = (zHigh-z)/(zHigh-zLow);
                        weightHigh = (z-zLow)/(zHigh-zLow);
                        for k=1:length(projIndexLow)
                            sgmHelical(j,:,i) = sgmHelical(j,:,i)+weightLow.*sgm(projIndexLow(k),:)+weightHigh.*sgm(projIndexHigh(k),:);
                        end
                    elseif (zLow>10.2346 && zLow<15.4734)
                        projIndexLow = find(zProjs>10.2346 & zProjs<10.4775); %Bass Ackwards way
                        projIndexHigh = find(zProjs>15.4734 & zProjs<15.7162);
                        weightLow = (zHigh-z)/(zHigh-zLow);
                        weightHigh = (z-zLow)/(zHigh-zLow);
                        for k=1:length(projIndexLow)
                            sgmHelical(j,:,i) = sgmHelical(j,:,i)+weightLow.*sgm(projIndexLow(k),:)+weightHigh.*sgm(projIndexHigh(k),:);
                        end
                    else
                        disp('Error: unused projection(s) at j==1 \n')
                    end
%%%
                end
           
            elseif j==scanFrames
                if (2.3765<z && z<2.6194) %inside rest zone and at angle 0
                    projIndexLow = find(zProjs>2.3765 & zProjs<2.6194); %Bass Ackwards way
                    for k=1:length(projIndexLow)
                        sgmHelical(j,:,i) = sgmHelical(j,:,i)+sgm(projIndexLow(k),:);
                    end
                elseif (7.6152<z && z<7.8581)
                    projIndexLow = find(zProjs>7.6152 & zProjs<7.8581); %Bass Ackwards way
                    for k=1:length(projIndexLow)
                        sgmHelical(j,:,i) = sgmHelical(j,:,i)+sgm(projIndexLow(k),:);
                    end
                elseif (12.8540<z && z<13.0969)
                    projIndexLow = find(zProjs>12.854 & zProjs<13.0969); %Bass Ackwards way
                    for k=1:length(projIndexLow)
                        sgmHelical(j,:,i) = sgmHelical(j,:,i)+sgm(projIndexLow(k),:);
                    end
                elseif (18.0927<z && z<18.3356) 
                    projIndexLow = find(zProjs>18.0927 & zProjs<18.3356); %Bass Ackwards way
                    for k=1:length(projIndexLow)
                        sgmHelical(j,:,i) = sgmHelical(j,:,i)+sgm(projIndexLow(k),:);
                    end
                else
                    projIndexLow = frameIndex+(sweepNumber)*(scanFrames+restFrames); % initial beta
                    zLow = zProjs(projIndexLow);
                    %360LI
                    projIndexHigh = projIndexLow+restFrames+2*(scanFrames-frameIndex)+1;
                    zHigh = zProjs(projIndexHigh);
                    
                    if zHigh < z % brute force and ignorance method of making sure its interpolation and not extrapolation
                        sweepNumber=sweepNumber+1;
                        frameIndex=scanFrames-frameIndex+1;
                        projIndexLow = frameIndex+(sweepNumber)*(scanFrames+restFrames); % initial beta
                        zLow = zProjs(projIndexLow);
                        projIndexHigh = projIndexLow+restFrames+2*(scanFrames-frameIndex)+1;
                        zHigh = zProjs(projIndexHigh);
                    end
%%%
                    if (zLow>2.3765 && zLow<7.6152)
                        projIndexLow = find(zProjs>2.3765 & zProjs<2.6194); %Bass Ackwards way
                        projIndexHigh = find(zProjs>7.6152 & zProjs<7.8581);
                        weightLow = (zHigh-z)/(zHigh-zLow);
                        weightHigh = (z-zLow)/(zHigh-zLow);
                        for k=1:length(projIndexLow)
                            sgmHelical(j,:,i) = sgmHelical(j,:,i)+weightLow.*sgm(projIndexLow(k),:)+weightHigh.*sgm(projIndexHigh(k),:);
                        end
                    elseif (zLow>7.6152 && zLow<12.8540)
                        projIndexLow = find(zProjs>7.6152 & zProjs<7.8581); %Bass Ackwards way
                        projIndexHigh = find(zProjs>12.8540 & zProjs<13.0969);
                        weightLow = (zHigh-z)/(zHigh-zLow);
                        weightHigh = (z-zLow)/(zHigh-zLow);
                        for k=1:length(projIndexLow)
                            sgmHelical(j,:,i) = sgmHelical(j,:,i)+weightLow.*sgm(projIndexLow(k),:)+weightHigh.*sgm(projIndexHigh(k),:);
                        end
                    elseif (zLow>12.8540 && zLow<18.0927)
                        projIndexLow = find(zProjs>12.8540 & zProjs<13.0969); %Bass Ackwards way
                        projIndexHigh = find(zProjs>18.0927 & zProjs<18.3356);
                        weightLow = (zHigh-z)/(zHigh-zLow);
                        weightHigh = (z-zLow)/(zHigh-zLow);
                        for k=1:length(projIndexLow)
                            sgmHelical(j,:,i) = sgmHelical(j,:,i)+weightLow.*sgm(projIndexLow(k),:)+weightHigh.*sgm(projIndexHigh(k),:);
                        end
                    else
                        disp('Error: unused projection(s) at j==1 \n')
                    end
%%%
                    
                end
            
            else % for all other projection angles
                projIndexLow = frameIndex+(sweepNumber)*(scanFrames+restFrames); % initial beta
                zLow = zProjs(projIndexLow);

                %360LI
                projIndexHigh = projIndexLow+restFrames+2*(scanFrames-frameIndex)+1;
                zHigh = zProjs(projIndexHigh);
                
                if zHigh < z % brute force and ignorance method of making sure its interpolation and not extrapolation
                    sweepNumber=sweepNumber+1;
                    frameIndex=scanFrames-frameIndex+1;
                    projIndexLow = frameIndex+(sweepNumber)*(scanFrames+restFrames); % initial beta
                    zLow = zProjs(projIndexLow);
                    projIndexHigh = projIndexLow+restFrames+2*(scanFrames-frameIndex)+1;
                    zHigh = zProjs(projIndexHigh);
                end
                
                weightLow = (zHigh-z)/(zHigh-zLow);
                weightHigh = (z-zLow)/(zHigh-zLow);
                sgmHelical(j,:,i) = weightLow.*sgm(projIndexLow,:)+weightHigh.*sgm(projIndexHigh,:);
            end

        end
    end
end

%%
filename = sprintf('sgm/20220405_ReverseHelical/360LIAVG_sgm_1_%d-%d-%d.raw', size(sgmHelical, 2), size(sgmHelical, 1), size(sgmHelical,3));
MgSaveRawFile(filename, sgmHelical);

%%
!mgfbp config_mgfbp_Treb_I_ReverseHelical.jsonc

% convert to mu
img = MgReadRawFile('img/ReverseHelical_tube_on/multislice/multislice_CIRS_120kV_2mA_img_bh_1_512-512-289.raw', 512, 512, 289);

RhoPMMA = 1.18; %g/cc
muOverRhoPMMA = 0.02067/RhoPMMA;
imgMU = img.*muOverRhoPMMA.*RhoPMMA; %convert to mu

%convert to HU
muWater = 0.019216;
imgHU = (imgMU-muWater)./muWater*1000;
MgSaveRawFile('img/ReverseHelical_tube_on/multislice/multislice_CIRS_HU_120kV_2mA_img_bh_1_512-512-289.raw', imgHU);
%% close all; clear; clc;
img = MgReadRawFile('img/ReverseHelical_tube_on/multislice/multislice_CIRS_120kV_2mA_img_bh_1_512-512-289.raw', 512, 512, 289);
for i=1:289
    imgtemp = img(:,:,i);
    %imgtemp = MgRingCorrection(imgtemp, -0.8, 100, 0.05, 30, 60);
    imgtemp = MgRingCorrection(imgtemp, -1000, 1000, 200, 30, 60);
    %imgtemp = MgRingCorrection(imgtemp, -200, 200, 200, 30, 70);
    img(:,:,i) = imgtemp;
end

MgSaveRawFile('img/ReverseHelical_tube_on/multislice/multislice_CIRS_120kV_2mA_img_bh_ring_1_512-512-289.raw', img);