clear
% 12/17/21

imgH = MgReadRawFile('img_1_HE.raw', 512,512,1);
imgL = MgReadRawFile('img_1_LE.raw', 512,512,1);

%convert HU images into mu images (only if input images are in HU)
imgH = (imgH./1000.*0.0148)+0.0148;
imgL = (imgL./1000.*0.015808)+0.015808;

%compute matrix
% MAT = [muWH, (muIH-muWH); muWL, (muIL-muWL)];
MAT = [0.018, 0.00361; 0.0178, 0.00596]; % matrix given by Mang
%MAT = [0.0156, 0.0023; 0.0165, 0.0033]; % My matrix

% apply inverse matrix to high and low images pixel by pixel
dim = length(imgH);
imgW = zeros(dim, dim);
imgI = imgW;
for i=1:dim
    for j=1:dim
         temp = MAT\[imgH(i,j); imgL(i,j)];
         imgW(i,j) = temp(1);
         imgI(i,j) = temp(2);
    end
end

MgSaveRawFile('img_water_iodine.raw', imgW)
MgSaveRawFile('img_iodine_water.raw', imgI)

%% correlation based denoising
imgW_FT = fftshift(fft2(imgW));
imgI_FT = fftshift(fft2(imgI));

[M N]=size(imgW); % image size
R=500; % filter size parameter 500
X=0:N-1;
Y=0:M-1;
[X Y]=meshgrid(X,Y);
Cx=0.5*N;
Cy=0.5*M;
Lo=exp(-((X-Cx).^2+(Y-Cy).^2)./(2*R).^2);
Hi=1-Lo; % High pass filter=1-low pass filter

imgW_HP = ifft2(ifftshift(imgW_FT.*Hi));
imgI_HP = ifft2(ifftshift(imgI_FT.*Hi));

MgSaveRawFile('img_water_iodine_HP.raw', imgW_HP)
MgSaveRawFile('img_iodine_water_HP.raw', imgI_HP)

%% trial and error weighting factors
for k=20 %empiricle weighting factor 20 for Mang matrix, 9 for my matrix
    denoiseW = imgW+k.*imgI_HP;
    %figure, imshow(denoiseW,[0.5, 1.1])
end

for k=190 %empiricle weighting factor 190 for Mang Matrix, 120 for my matrix
    denoiseI = imgI+k.*imgW_HP;
    %figure, imshow(denoiseI,[-0.1, 1.3])
end

MgSaveRawFile('img_water_iodine_DENOISE.raw', denoiseW)
MgSaveRawFile('img_iodine_water_DENOISE.raw', denoiseI)