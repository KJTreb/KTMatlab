close all; clear; clc;
%This function calls KTThorBadPixelCorr, KTThorPanelInterp, MgReadEviData,
%MgPolyval, MgWriteTiff; requires c1_PMMA_cali

%Fully processes 2D DE images from Thor PCD  

%load('20210824/coefficients/cali_PMMA.mat')
frames=1000;
%% TE

img_TE_data = MgReadEviData('20210824/Meat1_TE.EVI');
img_TE = img_TE_data(:,:,1:frames);
img_TE = mean(img_TE, 3);
clear img_TE_data;

air_TE = MgReadEviData('20210824/Calibration/PMMA=0_TE.EVI');
air_TE = mean(air_TE, 3);

% % bin pixels
% for i=1:255
%     for j=1:511
%         temp = img_TE(2*i:2*i+1,2*j:2*j+1);
%         img_TE_bin(i,j) = mean(temp(:));
%         temp = air_TE(2*i:2*i+1,2*j:2*j+1);
%         air_TE_bin(i,j) = mean(temp(:));
%     end
% end
% sino_TE = log(air_TE_bin./img_TE_bin);

% no pixel binning
sino_TE = log(air_TE./img_TE);

% %panel correction
% sgm_TE_corr = MgPolyval(cali_PMMA_TE, sino_TE);
% sgm_TE_corr(isnan(sgm_TE_corr)) = 0;
% sgm_TE_corr(isinf(sgm_TE_corr)) = 0;

%% HE

img_HE_data = MgReadEviData('20210824/Meat1_HE.EVI');
img_HE = img_HE_data(:,:,1:frames);
img_HE = mean(img_HE, 3);
clear img_HE_data;

air_HE = MgReadEviData('20210824/Calibration/PMMA=0_HE.EVI');
air_HE = mean(air_HE, 3);


% % bin pixels
% for i=1:255
%     for j=1:511
%         temp = img_HE(2*i:2*i+1,2*j:2*j+1);
%         img_HE_bin(i,j) = mean(temp(:));
%         temp = air_HE(2*i:2*i+1,2*j:2*j+1);
%         air_HE_bin(i,j) = mean(temp(:));
%     end
% end
% sino_HE = log(air_HE_bin./img_HE_bin);

% no pixel binning
sino_HE = log(air_HE./img_HE);

% % panel correction
% sgm_HE_corr = MgPolyval(cali_PMMA_HE, sino_HE);
% sgm_HE_corr(isnan(sgm_HE_corr)) = 0;
% sgm_HE_corr(isinf(sgm_HE_corr)) = 0;

%% LE
% % pixel binning
% img_LE_bin = img_TE_bin - img_HE_bin;
% air_LE_bin = air_TE_bin - air_HE_bin;
% sino_LE = log(air_LE_bin ./ img_LE_bin);

% no pixel binning
img_LE = img_TE - img_HE;
air_LE = air_TE - air_HE;
sino_LE = log(air_LE ./ img_LE);

% % panel correction
% sgm_LE_corr = MgPolyval(cali_PMMA_LE, sino_TE);
% sgm_LE_corr(isnan(sgm_LE_corr)) = 0;
% sgm_LE_corr(isinf(sgm_LE_corr)) = 0;

%% Write images
% %Interpolate known bad pixels
% sgm_TE_corr = KTThorBadPixelCorr(sgm_TE_corr, 1);
% sgm_HE_corr = KTThorBadPixelCorr(sgm_HE_corr, 1);
% sgm_LE_corr = KTThorBadPixelCorr(sgm_LE_corr, 1);

sino_TE = KTThorBadPixelCorr(sino_TE, 1);
sino_HE = KTThorBadPixelCorr(sino_HE, 1);
sino_LE = KTThorBadPixelCorr(sino_LE, 1);

% %Interpolate panel gaps
% sgm_TE_corr = KTThorPanelInterp(sgm_TE_corr);
% sgm_HE_corr = KTThorPanelInterp(sgm_HE_corr);
% sgm_LE_corr = KTThorPanelInterp(sgm_LE_corr);

sino_TE = KTThorPanelInterp(sino_TE);
sino_HE = KTThorPanelInterp(sino_HE);
sino_LE = KTThorPanelInterp(sino_LE);


MgWriteTiff('20210824/Processed/sgm_TE_no_corr.tif', sino_TE);
MgWriteTiff('20210824/Processed/sgm_HE_no_corr.tif', sino_HE);
MgWriteTiff('20210824/Processed/sgm_LE_no_corr.tif', sino_LE);
% MgWriteTiff('20210824/Processed/sgm_TE_corr.tif', sgm_TE_corr);
% MgWriteTiff('20210824/Processed/sgm_LE_corr.tif', sgm_LE_corr);
% MgWriteTiff('20210824/Processed/sgm_HE_corr.tif', sgm_HE_corr);

%% DE processing find weighting factor

% %20_85_meat (back=meat, obj(to be subtracted)=fat)
% HE_back = 0.44097;
% LE_back = 0.62135;
% HE_obj = 0.26662;
% LE_obj = 0.38133;
% 
% w_guess = (LE_obj-LE_back)/(HE_obj-HE_back);
% 
% % for no panel correction
% sino_TE(isnan(sino_TE)) = 0;
% sino_TE(isinf(sino_TE)) = 0;
% 
% sino_HE(isnan(sino_HE)) = 0;
% sino_HE(isinf(sino_HE)) = 0;
% 
% sino_LE(isnan(sino_LE)) = 0;
% sino_LE(isinf(sino_LE)) = 0;
% 
% for i=0.92 %for Gd subtraction
%    w = w_guess*i;
%    img_final = sgm_LE_corr-w.*sgm_HE_corr;
% %    img_final = sino_LE-w.*sino_HE;
%    figure, imshow(img_final,(1./i).*[mean(img_final(:))-3.*var(img_final(:)), mean(img_final(:))+3.*var(img_final(:))])
% end
% MgWriteTiff('20210820/20_85/Processed/WEIGHTED_SUBTRACTION_Gd.tif', img_final);
% 
% for i=1 %for Fat subtraction
%    w = w_guess*i;
%    img_final = sgm_LE_corr-w.*sgm_HE_corr;
% %    img_final = sino_LE-w.*sino_HE;
%    figure, imshow(img_final,(1./i).*[mean(img_final(:))-3.*var(img_final(:)), mean(img_final(:))+3.*var(img_final(:))])
% end
% MgWriteTiff('20210820/20_85/Processed/WEIGHTED_SUBTRACTION_fat.tif', img_final);
