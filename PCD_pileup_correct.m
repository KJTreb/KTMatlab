close all; clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate LUTs
res = 0.02; % LUT resolution: 0.02-> one count for 52.87 PCD slope

res_counts = 10; % PCD_true counts resolution in LUT (3 is good)
omega = 0.04; % percent of mA_fit values passed to EID LUT after selected for consistency with PCD counts (0.04 is "good")

%mA_fit = [res:res_counts*res:150]; %"fine" LUT
mA_fit = [0:res_counts*res:50, (50+5*res_counts*res):5*res_counts*res:400]; %"test" LUT

N=numel(mA_fit);
rows = 64;
cols = 5120;
LUT_PCD = zeros(rows,cols,N);


%raw EID signal
%LUT_EID = 7952.1./4.*mA_fit; % dividing by 4 to match signal/pixel for PCD pixel dimensions
% LUT_EID = (3489.605*30/120/16)./(0.5).*mA_fit; %Fluoro mode
LUT_EID = 3106/0.5.*mA_fit; %Fluoro mode

load('coefficients/p_PCD_HE_polynomial.mat');
for i=1:rows
    for j=1:cols
%         a=p_PCD_TE{i,j}(1);
%         b=p_PCD_TE{i,j}(2);
%         w=p_PCD_TE{i,j}(3);
%         LUT_PCD(i,j,:) = w.*(a.*mA_fit.*exp(-b.*a.*mA_fit))+(1-w).*(a.*mA_fit/(1+b.*a.*mA_fit));

        p1=p_PCD_HE{i,j}(1);
        p2=p_PCD_HE{i,j}(2);
        p3=p_PCD_HE{i,j}(3);
        p4=p_PCD_HE{i,j}(4);
        p5=p_PCD_HE{i,j}(5);
        p6=p_PCD_HE{i,j}(6);
        LUT_PCD(i,j,:) = p1.*mA_fit.^5 + p2.*mA_fit.^4 + p3.*mA_fit.^3 + p4.*mA_fit.^2 + p5.*mA_fit + p6;
    end
end


%LUT_PCD_true = 52.87.*mA_fit; %TE
LUT_PCD_true = 22.81.*mA_fit; %HE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
views = 1;
rows = 64;
cols = 5120;

EID_left = 495;
EID_right = 778;
EID_registered_width = 5120-EID_left-EID_right;

mA=10;
%t_PMMA = 9;
%prj_PCD = MgReadRawFile('PCDCT/data/air(old)/0.5mA_HE.raw',rows,cols,views);
prj_PCD = MgReadEviData('2D_PCD/data/star/10mA_TE.evi');


prj_EID_raw = MgReadTiff('2D_FPD_registered/data/star.tif');
prj_EID_scaled = prj_EID_raw.*mA./0.5; %compensate for fps, desired mA, and pixel size (fluoro mode)
for n = 1:views
    prj_EID(:,:,n) = [3401/0.5*mA.*ones(64, EID_left), prj_EID_scaled(:,:,n), 3196/0.5*mA.*ones(64, EID_right)]; %fill in area outside FPD FOV with air from LUT
end

% prj_EID = MgReadTiff('FPDCT_registered/Ca_final_new.tif');
% bkgnd = prj_EID(42, 42, 5);
% for n = 1:views
%     prj_EID(:,:,n) = [prj_EID(:,8:5120,n), bkgnd.*ones(64, 7)]; %fill in area outside FPD FOV with air from LUT
% end
% MgSaveRawFile('FPDCT_registered/test.raw', prj_EID);

%%% simulated data for now, later upsample registered FPD data
prj_EID = poissrnd(3106/0.5.*mA, rows,cols,views).*exp(-t_PMMA*0.235*1.18);
%%%

% prj_EID = MgReadRawFile('2D/10mA_1.raw',rows,cols,views,0,0,'uint16');
% prj_EID = prj_EID.*(mA/10).*(62120/600);
% prj_EID = imgaussfilt(prj_EID,1.5); %blur

M = numel(LUT_EID);
distance_PCD = zeros(1,M);
prj_PCD_corr = zeros(rows,cols,views);

% tic
% for k=1:views
%     for i=1:rows
%         for j=1:cols
%             N_PCD = prj_PCD(i,j,k);
%             N_EID = prj_EID(i,j,k);
%             for m=1:M
%                 distance(m) = ((LUT_PCD(i,j,m)-N_PCD)/N_PCD)^2 + omega.*((LUT_EID(m)-N_EID)/N_EID)^2; %L2 norm squared
%             end
%             [throwaway, idx_min] = max(-distance);
% %             if LUT_PCD_true(idx_min)<300 % don't do LUT when extremely low counts (reduces noise)
% %                 prj_PCD_corr(i,j,k) = (N_PCD+LUT_PCD_true(idx_min))./2;
% %             else
% %                 prj_PCD_corr(i,j,k) = LUT_PCD_true(idx_min); % corrected PCD projections (raw counts domain)
% %             end
%              prj_PCD_corr(i,j,k) = LUT_PCD_true(idx_min); 
%         end
%     end
%     time_elapsed = toc/60;
%     time_remaining = (toc/k)*(views-k)/60;
%     clc
%     disp(['Time elapsed: ',num2str(time_elapsed),' minutes. Time remaining: ',num2str(time_remaining),' minutes.'])
%     disp([num2str(k), '/', num2str(views)])
% end

tic
for k=1:views
    for i=1:rows
        for j=1:cols
            N_PCD = prj_PCD(i,j,k);
            N_EID = prj_EID(i,j,k);
            for m=1:M
                distance_PCD(m) = (LUT_PCD(i,j,m)-N_PCD)^2;
            end
            [throwaway, idx_mins] = mink(distance_PCD, ceil(omega*length(mA_fit))); % 4% of mA_fit values 
            distance = (LUT_EID(idx_mins)-N_EID).^2;
            [throwaway, idx_min] = max(-distance);

            %prj_PCD_corr(i,j,k) = LUT_PCD_true(idx_mins); 
            prj_PCD_corr(i,j,k) = LUT_PCD_true(idx_mins(idx_min)); 
        end
    end
    time_elapsed = toc/60;
    time_remaining = (toc/k)*(views-k)/60;
    clc
    disp(['Time elapsed: ',num2str(time_elapsed),' minutes. Time remaining: ',num2str(time_remaining),' minutes.'])
    disp([num2str(k), '/', num2str(views)])
end

MgSaveRawFile('PCDCT/PMMA/data_corr/PMMA=9_HE.raw', prj_PCD_corr);