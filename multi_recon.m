% 12/21/2021
clear
% NOTE: move images to a new location every time after recon so they dont get
% overwritten

% initialize variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EnergyBin = 'TE';
SliceCount = 10;
PixelBinning = 4;
% Also modify other variables in jsonc files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import data
if strcmp(EnergyBin, 'LE')
    muWater = 0.015808;
    objTE = MgReadEviDataCrop('data\1_TE.EVI'); % big long 10 sweep evi file
    objHE = MgReadEviDataCrop('data\1_HE.EVI'); % big long 10 sweep evi file
    obj = objTE-objHE;
    clear objTE
    clear objHE

    airTE = MgReadEviDataCrop('air\1_TE.EVI'); % air file for only one sweep
    airHE = MgReadEviDataCrop('air\1_HE.EVI'); % air file for only one sweep
    air = airTE-airHE;
    clear airTE
    clear airHE
elseif strcmp(EnergyBin, 'HE')
    muWater = 0.0148;
    objFilename = strcat('data\1_',EnergyBin,'.EVI');
    airFilename = strcat('air\1_',EnergyBin,'.EVI');
    obj = MgReadEviDataCrop(objFilename); % big long 10 sweep evi file
    air = MgReadEviDataCrop(airFilename); % air file for only one sweep
else 
    muWater = 0.015328; % TE
    objFilename = strcat('data\1_',EnergyBin,'.EVI');
    airFilename = strcat('air\1_',EnergyBin,'.EVI');
    obj = MgReadEviDataCrop(objFilename); % big long 10 sweep evi file
    air = MgReadEviDataCrop(airFilename); % air file for only one sweep
end
% generate sinogram
air=mean(air,3);
sgm=log(air./obj);
% sgm(isnan(sgm))=0;
% sgm(isinf(sgm))=0;
clear obj
clear air

% Split full sinogram into 10 sinograms (and flip reverse direction projections)
one=sgm(:,:,1:495);
two=sgm(:,:,496:2*495);
two=flip(two,3);
three=sgm(:,:,2*495+1:3*495);
four=sgm(:,:,3*495+1:4*495);
four=flip(four,3);
five=sgm(:,:,4*495+1:5*495);
six=sgm(:,:,5*495+1:6*495);
six=flip(six,3);
seven=sgm(:,:,6*495+1:7*495);
eight=sgm(:,:,7*495+1:8*495);
eight=flip(eight,3);
nine=sgm(:,:,8*495+1:9*495);
ten=sgm(:,:,9*495+1:10*495);
ten=flip(ten,3);
clear sgm

% set slice thickness
SliceThickness = 60/SliceCount;

one = MgResliceProjectionToSinogram(one, 3, SliceThickness, SliceCount);
two = MgResliceProjectionToSinogram(two, 3, SliceThickness, SliceCount);
three = MgResliceProjectionToSinogram(three, 3, SliceThickness, SliceCount);
four = MgResliceProjectionToSinogram(four, 3, SliceThickness, SliceCount);
five = MgResliceProjectionToSinogram(five, 3, SliceThickness, SliceCount);
six = MgResliceProjectionToSinogram(six, 3, SliceThickness, SliceCount);
seven = MgResliceProjectionToSinogram(seven, 3, SliceThickness, SliceCount);
eight = MgResliceProjectionToSinogram(eight, 3, SliceThickness, SliceCount);
nine = MgResliceProjectionToSinogram(nine, 3, SliceThickness, SliceCount);
ten = MgResliceProjectionToSinogram(ten, 3, SliceThickness, SliceCount);

%Panel correction 
one=permute(one, [2 1 3]);
one = XuPanelCorrection(one, 'PMMA', EnergyBin, SliceCount, 800, 4000, '7sDynamicCT-SinglePixel-th70-new');
one=permute(one, [2 1 3]);
two=permute(two, [2 1 3]);
two = XuPanelCorrection(two, 'PMMA', EnergyBin, SliceCount, 800, 4000, '7sDynamicCT-SinglePixel-th70-new');
two=permute(two, [2 1 3]);
three=permute(three, [2 1 3]);
three = XuPanelCorrection(three, 'PMMA', EnergyBin, SliceCount, 800, 4000, '7sDynamicCT-SinglePixel-th70-new');
three=permute(three, [2 1 3]);
four=permute(four, [2 1 3]);
four = XuPanelCorrection(four, 'PMMA', EnergyBin, SliceCount, 800, 4000, '7sDynamicCT-SinglePixel-th70-new');
four=permute(four, [2 1 3]);
five=permute(five, [2 1 3]);
five = XuPanelCorrection(five, 'PMMA', EnergyBin, SliceCount, 800, 4000, '7sDynamicCT-SinglePixel-th70-new');
five=permute(five, [2 1 3]);
six=permute(six, [2 1 3]);
six = XuPanelCorrection(six, 'PMMA', EnergyBin, SliceCount, 800, 4000, '7sDynamicCT-SinglePixel-th70-new');
six=permute(six, [2 1 3]);
seven=permute(seven, [2 1 3]);
seven = XuPanelCorrection(seven, 'PMMA', EnergyBin, SliceCount, 800, 4000, '7sDynamicCT-SinglePixel-th70-new');
seven=permute(seven, [2 1 3]);
eight=permute(eight, [2 1 3]);
eight = XuPanelCorrection(eight, 'PMMA', EnergyBin, SliceCount, 800, 4000, '7sDynamicCT-SinglePixel-th70-new');
eight=permute(eight, [2 1 3]);
nine=permute(nine, [2 1 3]);
nine = XuPanelCorrection(nine, 'PMMA', EnergyBin, SliceCount, 800, 4000, '7sDynamicCT-SinglePixel-th70-new');
nine=permute(nine, [2 1 3]);
ten=permute(ten, [2 1 3]);
ten = XuPanelCorrection(ten, 'PMMA', EnergyBin, SliceCount, 800, 4000, '7sDynamicCT-SinglePixel-th70-new');
ten=permute(ten, [2 1 3]);

% Pixel binning
sgmWidth = 5120/PixelBinning;

one=imresize(one,[495, sgmWidth],'box');
two=imresize(two,[495, sgmWidth],'box');
three=imresize(three,[495, sgmWidth],'box');
four=imresize(four,[495, sgmWidth],'box');
five=imresize(five,[495, sgmWidth],'box');
six=imresize(six,[495, sgmWidth],'box');
seven=imresize(seven,[495, sgmWidth],'box');
eight=imresize(eight,[495, sgmWidth],'box');
nine=imresize(nine,[495, sgmWidth],'box');
ten=imresize(ten,[495, sgmWidth],'box');

% truncation correction
pixelSize = PixelBinning/10;
left = (800/PixelBinning+1):(800/PixelBinning+16);
right = (4000/PixelBinning-16):(4000/PixelBinning-1);

one = MgTruncationSgmExtraploation(one,pixelSize, muWater, left, right);
two = MgTruncationSgmExtraploation(two,pixelSize, muWater, left, right);
three = MgTruncationSgmExtraploation(three,pixelSize, muWater, left, right);
four = MgTruncationSgmExtraploation(four,pixelSize, muWater, left, right);
five = MgTruncationSgmExtraploation(five,pixelSize,muWater, left, right);
six = MgTruncationSgmExtraploation(six,pixelSize, muWater, left, right);
seven = MgTruncationSgmExtraploation(seven,pixelSize,muWater, left, right);
eight = MgTruncationSgmExtraploation(eight,pixelSize, muWater, left, right);
nine = MgTruncationSgmExtraploation(nine,pixelSize, muWater, left, right);
ten = MgTruncationSgmExtraploation(ten,pixelSize, muWater, left, right);

% save files
MgSaveRawFile('sgm\one\sgm_1.raw',one);
MgSaveRawFile('sgm\two\sgm_1.raw',two);
MgSaveRawFile('sgm\three\sgm_1.raw',three);
MgSaveRawFile('sgm\four\sgm_1.raw',four);
MgSaveRawFile('sgm\five\sgm_1.raw',five);
MgSaveRawFile('sgm\six\sgm_1.raw',six);
MgSaveRawFile('sgm\seven\sgm_1.raw',seven);
MgSaveRawFile('sgm\eight\sgm_1.raw',eight);
MgSaveRawFile('sgm\nine\sgm_1.raw',nine);
MgSaveRawFile('sgm\ten\sgm_1.raw',ten);

%% recon one
jsFileName = 'config_fbp_one.jsonc';
XuModifyJsoncFile(jsFileName,'SliceCount', SliceCount);
XuModifyJsoncFile(jsFileName,'DetectorElementSize', pixelSize);
XuModifyJsoncFile(jsFileName,'SinogramWidth', sgmWidth);
!mgfbp config_fbp_one.jsonc

% ring correction after recon
js = MgReadJsoncFile(jsFileName);
filename = sprintf('./%s/img_1.raw', js.OutputDir);

img = MgReadRawFile(filename, js.ImageDimension, js.ImageDimension, SliceCount);
for i=1:js.SliceCount
    imgtemp = img(:,:,i);
    imgtemp = MgRingCorrection(imgtemp, -300, 630, 100, 20, 70);
    %circular crop image
    imgtemp = MgCropCircleFOV(imgtemp, js.ImageSize, -1000);
    img(:,:,i) = imgtemp;
end

filename = sprintf('./%s/img_1.raw', js.OutputDir);
MgSaveRawFile(filename, img);


%% recon two
jsFileName = 'config_fbp_two.jsonc';
XuModifyJsoncFile(jsFileName,'SliceCount', SliceCount);
XuModifyJsoncFile(jsFileName,'DetectorElementSize', pixelSize);
XuModifyJsoncFile(jsFileName,'SinogramWidth', sgmWidth);
!mgfbp config_fbp_two.jsonc

% ring correction after recon
js = MgReadJsoncFile(jsFileName);
filename = sprintf('./%s/img_1.raw', js.OutputDir);

img = MgReadRawFile(filename, js.ImageDimension, js.ImageDimension, SliceCount);
for i=1:js.SliceCount
    imgtemp = img(:,:,i);
    imgtemp = MgRingCorrection(imgtemp, -300, 630, 100, 20, 70);
    %circular crop image
    imgtemp = MgCropCircleFOV(imgtemp, js.ImageSize, -1000);
    img(:,:,i) = imgtemp;
end

filename = sprintf('./%s/img_1.raw', js.OutputDir);
MgSaveRawFile(filename, img);

%% recon three
jsFileName = 'config_fbp_three.jsonc';
XuModifyJsoncFile(jsFileName,'SliceCount', SliceCount);
XuModifyJsoncFile(jsFileName,'DetectorElementSize', pixelSize);
XuModifyJsoncFile(jsFileName,'SinogramWidth', sgmWidth);
!mgfbp config_fbp_three.jsonc

% ring correction after recon
js = MgReadJsoncFile(jsFileName);
filename = sprintf('./%s/img_1.raw', js.OutputDir);

img = MgReadRawFile(filename, js.ImageDimension, js.ImageDimension, SliceCount);
for i=1:js.SliceCount
    imgtemp = img(:,:,i);
    imgtemp = MgRingCorrection(imgtemp, -300, 630, 100, 20, 70);
    %circular crop image
    imgtemp = MgCropCircleFOV(imgtemp, js.ImageSize, -1000);
    img(:,:,i) = imgtemp;
end

filename = sprintf('./%s/img_1.raw', js.OutputDir);
MgSaveRawFile(filename, img);

%% recon four

jsFileName = 'config_fbp_four.jsonc';
XuModifyJsoncFile(jsFileName,'SliceCount', SliceCount);
XuModifyJsoncFile(jsFileName,'DetectorElementSize', pixelSize);
XuModifyJsoncFile(jsFileName,'SinogramWidth', sgmWidth);
!mgfbp config_fbp_four.jsonc

% ring correction after recon
js = MgReadJsoncFile(jsFileName);
filename = sprintf('./%s/img_1.raw', js.OutputDir);

img = MgReadRawFile(filename, js.ImageDimension, js.ImageDimension, SliceCount);
for i=1:js.SliceCount
    imgtemp = img(:,:,i);
    imgtemp = MgRingCorrection(imgtemp, -300, 630, 100, 20, 70);
    %circular crop image
    imgtemp = MgCropCircleFOV(imgtemp, js.ImageSize, -1000);
    img(:,:,i) = imgtemp;
end

filename = sprintf('./%s/img_1.raw', js.OutputDir);
MgSaveRawFile(filename, img);

%% recon five

jsFileName = 'config_fbp_five.jsonc';
XuModifyJsoncFile(jsFileName,'SliceCount', SliceCount);
XuModifyJsoncFile(jsFileName,'DetectorElementSize', pixelSize);
XuModifyJsoncFile(jsFileName,'SinogramWidth', sgmWidth);
!mgfbp config_fbp_five.jsonc

% ring correction after recon
js = MgReadJsoncFile(jsFileName);
filename = sprintf('./%s/img_1.raw', js.OutputDir);

img = MgReadRawFile(filename, js.ImageDimension, js.ImageDimension, SliceCount);
for i=1:js.SliceCount
    imgtemp = img(:,:,i);
    imgtemp = MgRingCorrection(imgtemp, -300, 630, 100, 20, 70);
    %circular crop image
    imgtemp = MgCropCircleFOV(imgtemp, js.ImageSize, -1000);
    img(:,:,i) = imgtemp;
end

filename = sprintf('./%s/img_1.raw', js.OutputDir);
MgSaveRawFile(filename, img);

%% recon six

jsFileName = 'config_fbp_six.jsonc';
XuModifyJsoncFile(jsFileName,'SliceCount', SliceCount);
XuModifyJsoncFile(jsFileName,'DetectorElementSize', pixelSize);
XuModifyJsoncFile(jsFileName,'SinogramWidth', sgmWidth);
!mgfbp config_fbp_six.jsonc

% ring correction after recon
js = MgReadJsoncFile(jsFileName);
filename = sprintf('./%s/img_1.raw', js.OutputDir);

img = MgReadRawFile(filename, js.ImageDimension, js.ImageDimension, SliceCount);
for i=1:js.SliceCount
    imgtemp = img(:,:,i);
    imgtemp = MgRingCorrection(imgtemp, -300, 630, 100, 20, 70);
    %circular crop image
    imgtemp = MgCropCircleFOV(imgtemp, js.ImageSize, -1000);
    img(:,:,i) = imgtemp;
end

filename = sprintf('./%s/img_1.raw', js.OutputDir);
MgSaveRawFile(filename, img);

%% recon seven

jsFileName = 'config_fbp_seven.jsonc';
XuModifyJsoncFile(jsFileName,'SliceCount', SliceCount);
XuModifyJsoncFile(jsFileName,'DetectorElementSize', pixelSize);
XuModifyJsoncFile(jsFileName,'SinogramWidth', sgmWidth);
!mgfbp config_fbp_seven.jsonc

% ring correction after recon
js = MgReadJsoncFile(jsFileName);
filename = sprintf('./%s/img_1.raw', js.OutputDir);

img = MgReadRawFile(filename, js.ImageDimension, js.ImageDimension, SliceCount);
for i=1:js.SliceCount
    imgtemp = img(:,:,i);
    imgtemp = MgRingCorrection(imgtemp, -300, 630, 100, 20, 70);
    %circular crop image
    imgtemp = MgCropCircleFOV(imgtemp, js.ImageSize, -1000);
    img(:,:,i) = imgtemp;
end

filename = sprintf('./%s/img_1.raw', js.OutputDir);
MgSaveRawFile(filename, img);

%% recon eight

jsFileName = 'config_fbp_eight.jsonc';
XuModifyJsoncFile(jsFileName,'SliceCount', SliceCount);
XuModifyJsoncFile(jsFileName,'DetectorElementSize', pixelSize);
XuModifyJsoncFile(jsFileName,'SinogramWidth', sgmWidth);
!mgfbp config_fbp_eight.jsonc

% ring correction after recon
js = MgReadJsoncFile(jsFileName);
filename = sprintf('./%s/img_1.raw', js.OutputDir);

img = MgReadRawFile(filename, js.ImageDimension, js.ImageDimension, SliceCount);
for i=1:js.SliceCount
    imgtemp = img(:,:,i);
    imgtemp = MgRingCorrection(imgtemp, -300, 630, 100, 20, 70);
    %circular crop image
    imgtemp = MgCropCircleFOV(imgtemp, js.ImageSize, -1000);
    img(:,:,i) = imgtemp;
end

filename = sprintf('./%s/img_1.raw', js.OutputDir);
MgSaveRawFile(filename, img);

%% recon nine

jsFileName = 'config_fbp_nine.jsonc';
XuModifyJsoncFile(jsFileName,'SliceCount', SliceCount);
XuModifyJsoncFile(jsFileName,'DetectorElementSize', pixelSize);
XuModifyJsoncFile(jsFileName,'SinogramWidth', sgmWidth);
!mgfbp config_fbp_nine.jsonc

% ring correction after recon
js = MgReadJsoncFile(jsFileName);
filename = sprintf('./%s/img_1.raw', js.OutputDir);

img = MgReadRawFile(filename, js.ImageDimension, js.ImageDimension, SliceCount);
for i=1:js.SliceCount
    imgtemp = img(:,:,i);
    imgtemp = MgRingCorrection(imgtemp, -300, 630, 100, 20, 70);
    %circular crop image
    imgtemp = MgCropCircleFOV(imgtemp, js.ImageSize, -1000);
    img(:,:,i) = imgtemp;
end

filename = sprintf('./%s/img_1.raw', js.OutputDir);
MgSaveRawFile(filename, img);

%% recon ten

jsFileName = 'config_fbp_ten.jsonc';
XuModifyJsoncFile(jsFileName,'SliceCount', SliceCount);
XuModifyJsoncFile(jsFileName,'DetectorElementSize', pixelSize);
XuModifyJsoncFile(jsFileName,'SinogramWidth', sgmWidth);
!mgfbp config_fbp_ten.jsonc

% ring correction after recon
js = MgReadJsoncFile(jsFileName);
filename = sprintf('./%s/img_1.raw', js.OutputDir);

img = MgReadRawFile(filename, js.ImageDimension, js.ImageDimension, SliceCount);
for i=1:js.SliceCount
    imgtemp = img(:,:,i);
    imgtemp = MgRingCorrection(imgtemp, -300, 630, 100, 20, 70);
    %circular crop image
    imgtemp = MgCropCircleFOV(imgtemp, js.ImageSize, -1000);
    img(:,:,i) = imgtemp;
end
filename = sprintf('./%s/img_1.raw', js.OutputDir);
MgSaveRawFile(filename, img);