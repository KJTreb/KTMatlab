close all; clear; clc;

air = MgReadEviData("data/20210813_Screw1/air_TE.evi");
obj = MgReadEviData("data/20210813_Screw1/1_TE.evi");

air=KTThorBadPixelCorr(air);
obj=KTThorBadPixelCorr(obj);

air=mean(air,3);

sgm=log(air./obj);

sgm(isnan(sgm))=0;
sgm(isinf(sgm))=0;

sgm = MgResliceProjectionToSinogram(sgm, 1, 1, 512); %512 slices, 0.07 mm isotropic with 1024X1024

%sgm = MgResliceProjectionToSinogram(sgm, 1, 2, 256); %256 slices, 0.14 mm isotropic with 512X512

%sgm = MgResliceProjectionToSinogram(sgm, 1, 512, 1); %1 slice

%sgm = MgResliceProjectionToSinogram(sgm, 2, 6, 85); % 85 slices, 0.47 mm axial and 0.42 z pixel sizes (match to other PCD Carm) with 153X153
%sgm = MgRebinSinogram(sgm,4);

%Truncation correction (change fbp config file before running this)
sgm = XuExtrapSgmCarm_KT_Benchtop(sgm,'config_pre_temp','config_mgfbp_cone_beam_full_scan.jsonc');

%sgm = XuThorPanelInterp(sgm);
filename = sprintf('sgm/20210813_Screw1/sgm_SP_1_%d-%d-%d.raw', size(sgm, 2), size(sgm, 1), size(sgm, 3));
MgSaveRawFile(filename, sgm);
%%
!mgfbp config_mgfbp_cone_beam_full_scan.jsonc
