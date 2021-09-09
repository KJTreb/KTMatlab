close all; clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read PMMA calibraiton data and take log
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slabs = [0, 3, 6, 9, 12, 15];
N = numel(slabs);
rows = 512;
cols = 1024;
prj_TE = zeros(rows, cols, N);
prj_HE = zeros(rows, cols, N);

for k = 1:N
    filename = sprintf('20210820/20_85/Calibration/PMMA=%d_TE.EVI', slabs(k));
    temp = MgReadEviData(filename);
    prj_TE(:,:,k) = mean(temp,3);
    clear temp
    filename = sprintf('20210820/20_85/Calibration/PMMA=%d_HE.EVI', slabs(k));
    temp = MgReadEviData(filename);
    prj_HE(:,:,k) = mean(temp,3);
    clear temp
end

prj_LE = prj_TE - prj_HE;
prj_LE(prj_LE<0) = 0;

% Take log
sgm_TE = log(prj_TE(:,:,1) ./ prj_TE);
sgm_TE(isnan(sgm_TE)) = 0;
sgm_TE(isinf(sgm_TE)) = 0;

sgm_LE = log(prj_LE(:,:,1) ./ prj_LE);
sgm_LE(isnan(sgm_LE)) = 0;
sgm_LE(isinf(sgm_LE)) = 0;

sgm_HE = log(prj_HE(:,:,1) ./ prj_HE);
sgm_HE(isnan(sgm_HE)) = 0;
sgm_HE(isinf(sgm_HE)) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the calibration calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off','all')

n = 3;

% calculate the position of each detector pixel
offcenter_axial = 7.5;
u = (1:cols) - (cols/2) - 0.5;
u = u*0.1 + offcenter_axial;

offcenter_z = 0;
v = (1:rows) - (rows/2) - 0.5;
v = v*0.1 + offcenter_z;


% PMMA thickness
thickness = slabs / 2; % unit: inch
thickness = thickness * 25.4;

% system geometry
sdd = 1041;

cali_PMMA_TE = cell(rows, cols);
cali_PMMA_HE = cell(rows, cols);
cali_PMMA_LE = cell(rows, cols);

pb = MgCmdLineProgressBar('Calibrating: ');
for col = 1:cols
    pb.print(col, cols);

    for row = 1:rows
        % calculate pass length
        pass_length = zeros(1, N);
        for t = 1:N
            pass_length(t) = thickness(t)/sdd * sqrt(sdd^2 + u(col)^2 +v(row)^2);
        end
        % get projection data
        % LE
        p_LE = sgm_LE(row, col, :);
        p_LE = squeeze(p_LE)';
        % HE
        p_HE = sgm_HE(row, col, :);
        p_HE = squeeze(p_HE)';
        % TE
        p_TE = sgm_TE(row, col, :);
        p_TE = squeeze(p_TE)';
        
        % polynomial fit
        cali_PMMA_TE{row, col} = polyfit(p_TE, pass_length, n);
        cali_PMMA_LE{row, col} = polyfit(p_LE, pass_length, n);
        cali_PMMA_HE{row, col} = polyfit(p_HE(1:end-1), pass_length(1:end-1), n);
        
    end
end

% save parameters to file
MgMkdir('20210820/20_85/coefficients', false);
save('20210820/20_85/coefficients/cali_PMMA', 'cali_PMMA_HE', 'cali_PMMA_LE', 'cali_PMMA_TE');








