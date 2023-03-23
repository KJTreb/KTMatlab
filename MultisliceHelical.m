close all; clear; clc;

%%%%%%%%%%360LI algorithm for reverse extended helical trajectory, benchtop, PCD operated as if it only had a
%%%%%%%%%%single row

%Initialize parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
detz = 4.3; %detector z-coverage
%pitch = 0.6; % for shortscan plus rest
fps = 10;
%omega = 4.16; %deg/s
totFrames = 4350;
scanFrames = 500; %frames per short scan
restFrames = 50; %frames per rest period
numSweeps = 8; %scan range = numSweeps*v
slices = 60; %recon slices per detz, 14 to match FPD with 2x2 binning
PixelBinning = 4;
v3 = 0.00375;
zPerFrame = v3*25.4/fps; % mm/frame

zScan = scanFrames*zPerFrame; %obj z-translation (mm) per short scan
zRest = restFrames*zPerFrame; %obj z-translation (mm) per rest period

%extra 180 LI parameters
% DetectorOffcenter = 43;
% SDD = 1041.4;

%%
%Import data and generate initial sinogram with corrections %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
objFilename = 'data\20221114_Defrise\3_TE.EVI';
airFilename = 'data\20221114_Defrise\air_TE.EVI';

%obj = MgReadEviData(objFilename);
obj = MgReadRawFile('data\20221114_Defrise\4_TE.raw', 64, 5120, 4400, 0, 0, 'uint16');
air = MgReadEviData(airFilename);

air = mean(air, 3);

temp = air./obj;
clear air
clear obj
clear airFilename
clear objFilename

sgm = log(temp);
clear temp

sgm(isnan(sgm)) = 0;
sgm(isinf(sgm)) = 0;

%%
% %Panel correction
% load('coefficients120kV2mA.mat');
% sgm = MgPolyval(cali_PMMA_TE, sgm);
% sgm(isnan(sgm)) = 0;
% sgm(isinf(sgm)) = 0;
% %%%%%
% 
% %Panel gap interpolation
% sgm = MgInterpSinogram(sgm, 2, 2, 4);

sgm = MgResliceProjectionToSinogram(sgm, 3, 1, 60);%I module PCD 60 slices
sgm = MgRebinSinogram(sgm,PixelBinning); %pixel binning

sgmOG=sgm;
%% 4/4 360LI reverse helical trajectory
% notes, extrap z conditions for if statements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%known delay between data acquisition and start of phantom rotation
delay=10; %(36) can be up to 50, if higher then reaquire
sgmStageTwo=sgmOG((delay+1):(4400-(50-delay)),:,:);

% empirical offset between start of rotation and start of fram acquisition
% (up to 0.5 times frame time)
sgm = zeros(size(sgmStageTwo));
offset = 0.0; % XXX/fps (empirical, between 0 and 1, used as weighting factor)
for p=1:(length(sgm(:,1))-1)
    sgm(p,:,:) = (1-offset)*sgmStageTwo(p,:,:) + offset*sgmStageTwo(p+1,:,:);
end
%note don't need to account for reverse sweep direction because we are
%working with the original full-acquisition sinogram prior to "flipping"
%the direction in the construction of the final sinogram used for recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sgmHelical=zeros(scanFrames,5120/PixelBinning,floor(zPerFrame*totFrames*slices/detz));
zProjs=0:zPerFrame:zPerFrame*(totFrames-1);

% thetaProjs1 = 0:(208/scanFrames):(208/scanFrames*(scanFrames-1));
% thetaProjs2 = [thetaProjs1,208/scanFrames*(scanFrames-1)*ones(1,50)];
% thetaProjs3 = (208/scanFrames*(scanFrames-1)):(-208/scanFrames):0;
% thetaProjs4 = cat(2, thetaProjs2, thetaProjs3, zeros(1,50));
% 
% thetaProjs = cat(2, thetaProjs4, thetaProjs4, thetaProjs4, thetaProjs2, thetaProjs3); %in degrees, source angle between 0 and ~208

for i=1:(floor(zPerFrame*totFrames*slices/detz)-1) % through recon z
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
            
%             if j==1
%                 if (4.9959<z && z<5.2387) %inside rest zone and at angle 0
%                     projIndexLow = find(zProjs>4.9959 & zProjs<5.2387); %Bass Ackwards way
%                     for k=1:length(projIndexLow)
%                         sgmHelical(j,:,i) = sgmHelical(j,:,i)+sgm(projIndexLow(k),:);
%                     end
%                 elseif (10.2346<z && z<10.4775)
%                     projIndexLow = find(zProjs>10.2346 & zProjs<10.4775); %Bass Ackwards way
%                     for k=1:length(projIndexLow)
%                         sgmHelical(j,:,i) = sgmHelical(j,:,i)+sgm(projIndexLow(k),:);
%                     end
%                 elseif (15.4734<z && z<15.7162)
%                     projIndexLow = find(zProjs>15.4734 & zProjs<15.7162); %Bass Ackwards way
%                     for k=1:length(projIndexLow)
%                         sgmHelical(j,:,i) = sgmHelical(j,:,i)+sgm(projIndexLow(k),:);
%                     end
%                 else
%                     projIndexLow = frameIndex+(sweepNumber)*(scanFrames+restFrames); % initial beta
%                     zLow = zProjs(projIndexLow);
%                     %360LI
%                     projIndexHigh = projIndexLow+restFrames+2*(scanFrames-frameIndex)+1; % minus rest frames ???
%                     zHigh = zProjs(projIndexHigh);
%                     
%                     if zHigh < z % brute force and ignorance method of making sure its interpolation and not extrapolation
%                         sweepNumber=sweepNumber+1;
%                         frameIndex=scanFrames-frameIndex+1;
%                         projIndexLow = frameIndex+(sweepNumber)*(scanFrames+restFrames); % initial beta
%                         zLow = zProjs(projIndexLow);
%                         projIndexHigh = projIndexLow+restFrames+2*(scanFrames-frameIndex)+1;
%                         zHigh = zProjs(projIndexHigh);
%                     end
% %%%
%                     if zLow<4.9959
%                         projIndexHigh = find(zProjs>4.9959 & zProjs<5.2387);
%                         weightLow = (zHigh-z)/(zHigh-zLow);
%                         weightHigh = (z-zLow)/(zHigh-zLow);
%                         for k=1:length(projIndexHigh)
%                             sgmHelical(j,:,i) = sgmHelical(j,:,i)+weightLow.*sgm(projIndexLow,:)./length(projIndexHigh)+weightHigh.*sgm(projIndexHigh(k),:);
%                         end
%                     elseif (zLow>4.9959 && zLow<10.2346)
%                         projIndexLow = find(zProjs>4.9959 & zProjs<5.2387); %Bass Ackwards way
%                         projIndexHigh = find(zProjs>10.2346 & zProjs<10.4775);
%                         weightLow = (zHigh-z)/(zHigh-zLow);
%                         weightHigh = (z-zLow)/(zHigh-zLow);
%                         for k=1:length(projIndexLow)
%                             sgmHelical(j,:,i) = sgmHelical(j,:,i)+weightLow.*sgm(projIndexLow(k),:)+weightHigh.*sgm(projIndexHigh(k),:);
%                         end
%                     elseif (zLow>10.2346 && zLow<15.4734)
%                         projIndexLow = find(zProjs>10.2346 & zProjs<10.4775); %Bass Ackwards way
%                         projIndexHigh = find(zProjs>15.4734 & zProjs<15.7162);
%                         weightLow = (zHigh-z)/(zHigh-zLow);
%                         weightHigh = (z-zLow)/(zHigh-zLow);
%                         for k=1:50%length(projIndexLow)
%                             sgmHelical(j,:,i) = sgmHelical(j,:,i)+weightLow.*sgm(projIndexLow(k),:)+weightHigh.*sgm(projIndexHigh(k),:);
%                         end
%                     else
%                         disp('Error: unused projection(s) at j==1 \n')
%                     end
% %%%
% 
%                 end
%            
%             elseif j==scanFrames
%                 if (2.3765<z && z<2.6194) %inside rest zone and at angle 0
%                     projIndexLow = find(zProjs>2.3765 & zProjs<2.6194); %Bass Ackwards way
%                     for k=1:length(projIndexLow)
%                         sgmHelical(j,:,i) = sgmHelical(j,:,i)+sgm(projIndexLow(k),:);
%                     end
%                 elseif (7.6152<z && z<7.8581)
%                     projIndexLow = find(zProjs>7.6152 & zProjs<7.8581); %Bass Ackwards way
%                     for k=1:length(projIndexLow)
%                         sgmHelical(j,:,i) = sgmHelical(j,:,i)+sgm(projIndexLow(k),:);
%                     end
%                 elseif (12.8540<z && z<13.0969)
%                     projIndexLow = find(zProjs>12.854 & zProjs<13.0969); %Bass Ackwards way
%                     for k=1:length(projIndexLow)
%                         sgmHelical(j,:,i) = sgmHelical(j,:,i)+sgm(projIndexLow(k),:);
%                     end
%                 elseif (18.0927<z && z<18.3356) 
%                     projIndexLow = find(zProjs>18.0927 & zProjs<18.3356); %Bass Ackwards way
%                     for k=1:length(projIndexLow)
%                         sgmHelical(j,:,i) = sgmHelical(j,:,i)+sgm(projIndexLow(k),:);
%                     end
%                 else
%                     projIndexLow = frameIndex+(sweepNumber)*(scanFrames+restFrames); % initial beta
%                     zLow = zProjs(projIndexLow);
%                     %360LI
%                     projIndexHigh = projIndexLow+restFrames+2*(scanFrames-frameIndex)+1;
%                     zHigh = zProjs(projIndexHigh);
%                     
%                     if zHigh < z % brute force and ignorance method of making sure its interpolation and not extrapolation
%                         sweepNumber=sweepNumber+1;
%                         frameIndex=scanFrames-frameIndex+1;
%                         projIndexLow = frameIndex+(sweepNumber)*(scanFrames+restFrames); % initial beta
%                         zLow = zProjs(projIndexLow);
%                         projIndexHigh = projIndexLow+restFrames+2*(scanFrames-frameIndex)+1;
%                         zHigh = zProjs(projIndexHigh);
%                     end
%%%
%                     if (zLow>2.3765 && zLow<7.6152)
%                         projIndexLow = find(zProjs>2.3765 & zProjs<2.6194); %Bass Ackwards way
%                         projIndexHigh = find(zProjs>7.6152 & zProjs<7.8581);
%                         weightLow = (zHigh-z)/(zHigh-zLow);
%                         weightHigh = (z-zLow)/(zHigh-zLow);
%                         for k=1:length(projIndexLow)
%                             sgmHelical(j,:,i) = sgmHelical(j,:,i)+weightLow.*sgm(projIndexLow(k),:)+weightHigh.*sgm(projIndexHigh(k),:);
%                         end
%                     elseif (zLow>7.6152 && zLow<12.8540)
%                         projIndexLow = find(zProjs>7.6152 & zProjs<7.8581); %Bass Ackwards way
%                         projIndexHigh = find(zProjs>12.8540 & zProjs<13.0969);
%                         weightLow = (zHigh-z)/(zHigh-zLow);
%                         weightHigh = (z-zLow)/(zHigh-zLow);
%                         for k=1:length(projIndexLow)
%                             sgmHelical(j,:,i) = sgmHelical(j,:,i)+weightLow.*sgm(projIndexLow(k),:)+weightHigh.*sgm(projIndexHigh(k),:);
%                         end
%                     elseif (zLow>12.8540 && zLow<18.0927)
%                         projIndexLow = find(zProjs>12.8540 & zProjs<13.0969); %Bass Ackwards way
%                         projIndexHigh = find(zProjs>18.0927 & zProjs<18.3356);
%                         weightLow = (zHigh-z)/(zHigh-zLow);
%                         weightHigh = (z-zLow)/(zHigh-zLow);
%                         for k=1:length(projIndexLow)
%                             sgmHelical(j,:,i) = sgmHelical(j,:,i)+weightLow.*sgm(projIndexLow(k),:)+weightHigh.*sgm(projIndexHigh(k),:);
%                         end
%                     else
%                         disp('Error: unused projection(s) at j==scanFrames \n')
%                     end
%%%
                   
%                 end
            
            %else % for all other projection angles
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
                sliceLocLow = 30-(z-zLow)/0.1; 
                sliceLocHigh = 30+(zHigh-z)/0.1; 
                if 1<=sliceLocLow && sliceLocLow<=60 && 1<=sliceLocHigh && sliceLocHigh<=60
                    sgmHelical(j,:,i) = weightLow.*sgm(projIndexLow,:,round(sliceLocLow))+weightHigh.*sgm(projIndexHigh,:,round(sliceLocHigh)); 
                %%%%%%%%%%%%%%%%%%%%%%%%%
%                 elseif sliceLocLow<1 && sliceLocHigh>60
%                     sliceLocLow = 1;
%                     sliceLocHigh = 60;
                %%%%%%%%%%%%%%%%%%%%%%%%%
                elseif (z-zLow)<(zHigh-z)
                    sliceLoc = 30-(z-zLow)/0.1; 
                    if sliceLoc<0
                        sgmHelical(j,:,i) = sgm(projIndexLow,:,1);
                    else
                        sgmHelical(j,:,i) = sgm(projIndexLow,:,ceil(sliceLoc));
                    end
                else
                    sliceLoc = 30+(zHigh-z)/0.1; 
                    if sliceLoc>61
                        sgmHelical(j,:,i) = sgm(projIndexHigh,:,60);
                    else
                        sgmHelical(j,:,i) = sgm(projIndexHigh,:,floor(sliceLoc));
                    end
                end
            %end

        end
    end
end

%debugging
%sgmHelical(1,:,:) = 100000000*sgmHelical(1,:,:);
%sgmHelical = cat(1,sgmHelical,zeros(365,1280,67));

filename = sprintf('sgm/20221114_Defrise/sgm_4_%d-%d-%d.raw', size(sgmHelical, 2), size(sgmHelical, 1), size(sgmHelical,3));
MgSaveRawFile(filename, sgmHelical);

%%
!mgfbp config_mgfbp_Treb_I_ReverseHelical.jsonc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scale final image to account for upweighting of "rest time" projections
% img = MgReadRawFile('img/ReverseHelical_tube_on/head_lowmA_48point9_img_bh_1_512-512-67.raw', 512, 512, 67);
% 
% img = img.*500/550; %downweight final image
% 
% MgSaveRawFile('img/ReverseHelical_tube_on/head_lowmA_48point9_img_bh_1_512-512-67.raw', img);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% % close all; clear; clc;
% img = MgReadRawFile('img/20220324_QA_center_narrow/img_bh_1_512-512-60.raw', 512, 512, 60);
% for i=1:60
%     imgtemp = img(:,:,i);
%     imgtemp = MgRingCorrection(imgtemp, -1000000, 1000000, 100, 20, 70);
%     img(:,:,i) = imgtemp;
% end
% 
% MgSaveRawFile('img/20220324_QA_center_narrow/img_bh_ring_1_512-512-60.raw', img);