function [CRLB_I, CRLB_Gd, CRLB_PMMA] = Get_CRLB_four_bins(spectrum_Low, spectrum_High, PCD_threshold)
%% initialize spectra (already normalized)

%spectrum_Low = readmatrix('80kV.xlsx');
%spectrum_High = readmatrix('120kV.xlsx');

counts_High = spectrum_High(:,2);
counts_Low = spectrum_Low(:,2);

energy_High = spectrum_High(:,1);
energy_Low = spectrum_Low(:,1);
%% calculate parameters for CRLB calculations

temp = readmatrix('mu_I.xlsx');
mu_rho_1 = interp1(temp(:,1)*1000, temp(:,2), energy_High, 'linear'); %vector of mass attenuation coefficients for same energies as the spectra, units cm^2/gram
mu_1 = mu_rho_1;

temp = readmatrix('mu_Gd.xlsx');
mu_rho_2 = interp1(temp(:,1)*1000, temp(:,2), energy_High, 'linear'); %vector of mass attenuation coefficients for same energies as the spectra, units cm^2/gram
mu_2 = mu_rho_2;

temp = readmatrix('mu_PMMA.xlsx');
mu_rho_3 = interp1(temp(:,1)*1000, temp(:,2), energy_High, 'linear'); %vector of mass attenuation coefficients for same energies as the spectra, units cm^2/gram
mu_3 = mu_rho_3;


rho_1 = 0.010; % concentration in g/cc
rho_2 = 0.010; % concentration in g/cc
rho_3 = 1.18; % concentration in g/cc

A_1 = 2.5*rho_1; %iodine solution thickness in cm
A_2 = 2.5*rho_2; %Gd solution thickness in cm
A_3 = 16*rho_3; %PMMA thickness in cm

% get post-object spectra
% note, PMMA= 16 cm is already included in the imported spectra
counts_High = counts_High.*exp(-mu_1*A_1-mu_2*A_2);
counts_Low = [counts_Low; zeros(111-length(counts_Low),1)];
counts_Low = counts_Low.*exp(-mu_1*A_1-mu_2*A_2);



%% Optimize PCD threshold (note, assuming that mAs levels are constant for each spectrum)
index = 1;
for PCD_thresh = PCD_threshold
    
    % get PCD response spectra post-object
    bin_spectrum_High_LE = MgGetBinResponseSpectrum(counts_High, energy_High, 20, PCD_thresh);
    bin_spectrum_High_HE = MgGetBinResponseSpectrum(counts_High, energy_High, PCD_thresh, 120);
    bin_spectrum_Low_LE = MgGetBinResponseSpectrum(counts_Low, energy_High, 20, PCD_thresh);
    bin_spectrum_Low_HE = MgGetBinResponseSpectrum(counts_Low, energy_High, PCD_thresh, 120);
    
    % combine acquired data to get energy "bins" used for image formation
    bin_spectrum_1 = bin_spectrum_Low_LE;
    bin_spectrum_2 = bin_spectrum_High_LE;
    bin_spectrum_3 = bin_spectrum_Low_HE;
    bin_spectrum_4 = bin_spectrum_High_HE;
    
    %average number of photons detected for each PCD energy bin/spectrum combo
    %(expected or measured) through 16 cm PMMA
    N_High = 32.96; % 120 kV spectrum TE bin
    N_Low = 25.25; % low kV spectrum TE bin
    
    N_1=N_Low*sum(bin_spectrum_Low_LE);
    N_2=N_High*sum(bin_spectrum_High_LE);
    N_3=N_Low*sum(bin_spectrum_Low_HE);
    N_4=N_High*sum(bin_spectrum_High_HE);

    % Calculate Fisher information matrix
    M = zeros(3,3);
    M(1,1)=N_1.*sum(mu_1.*bin_spectrum_1)/sum(bin_spectrum_1).*sum(mu_1.*bin_spectrum_1)/sum(bin_spectrum_1)+N_2.*sum(mu_1.*bin_spectrum_2)/sum(bin_spectrum_2).*sum(mu_1.*bin_spectrum_2)/sum(bin_spectrum_2)+N_3.*sum(mu_1.*bin_spectrum_3)/sum(bin_spectrum_3).*sum(mu_1.*bin_spectrum_3)/sum(bin_spectrum_3)+N_4.*sum(mu_1.*bin_spectrum_4)/sum(bin_spectrum_4).*sum(mu_1.*bin_spectrum_4)/sum(bin_spectrum_4);
    M(1,2)=N_1.*sum(mu_1.*bin_spectrum_1)/sum(bin_spectrum_1).*sum(mu_2.*bin_spectrum_1)/sum(bin_spectrum_1)+N_2.*sum(mu_1.*bin_spectrum_2)/sum(bin_spectrum_2).*sum(mu_2.*bin_spectrum_2)/sum(bin_spectrum_2)+N_3.*sum(mu_1.*bin_spectrum_3)/sum(bin_spectrum_3).*sum(mu_2.*bin_spectrum_3)/sum(bin_spectrum_3)+N_4.*sum(mu_1.*bin_spectrum_4)/sum(bin_spectrum_4).*sum(mu_2.*bin_spectrum_4)/sum(bin_spectrum_4);
    M(1,3)=N_1.*sum(mu_1.*bin_spectrum_1)/sum(bin_spectrum_1).*sum(mu_3.*bin_spectrum_1)/sum(bin_spectrum_1)+N_2.*sum(mu_1.*bin_spectrum_2)/sum(bin_spectrum_2).*sum(mu_3.*bin_spectrum_2)/sum(bin_spectrum_2)+N_3.*sum(mu_1.*bin_spectrum_3)/sum(bin_spectrum_3).*sum(mu_3.*bin_spectrum_3)/sum(bin_spectrum_3)+N_4.*sum(mu_1.*bin_spectrum_4)/sum(bin_spectrum_4).*sum(mu_3.*bin_spectrum_4)/sum(bin_spectrum_4);
    M(2,1)=N_1.*sum(mu_2.*bin_spectrum_1)/sum(bin_spectrum_1).*sum(mu_1.*bin_spectrum_1)/sum(bin_spectrum_1)+N_2.*sum(mu_2.*bin_spectrum_2)/sum(bin_spectrum_2).*sum(mu_1.*bin_spectrum_2)/sum(bin_spectrum_2)+N_3.*sum(mu_2.*bin_spectrum_3)/sum(bin_spectrum_3).*sum(mu_1.*bin_spectrum_3)/sum(bin_spectrum_3)+N_4.*sum(mu_2.*bin_spectrum_4)/sum(bin_spectrum_4).*sum(mu_1.*bin_spectrum_4)/sum(bin_spectrum_4);
    M(2,2)=N_1.*sum(mu_2.*bin_spectrum_1)/sum(bin_spectrum_1).*sum(mu_2.*bin_spectrum_1)/sum(bin_spectrum_1)+N_2.*sum(mu_2.*bin_spectrum_2)/sum(bin_spectrum_2).*sum(mu_2.*bin_spectrum_2)/sum(bin_spectrum_2)+N_3.*sum(mu_2.*bin_spectrum_3)/sum(bin_spectrum_3).*sum(mu_2.*bin_spectrum_3)/sum(bin_spectrum_3)+N_4.*sum(mu_2.*bin_spectrum_4)/sum(bin_spectrum_4).*sum(mu_2.*bin_spectrum_4)/sum(bin_spectrum_4);
    M(2,3)=N_1.*sum(mu_2.*bin_spectrum_1)/sum(bin_spectrum_1).*sum(mu_3.*bin_spectrum_1)/sum(bin_spectrum_1)+N_2.*sum(mu_2.*bin_spectrum_2)/sum(bin_spectrum_2).*sum(mu_3.*bin_spectrum_2)/sum(bin_spectrum_2)+N_3.*sum(mu_2.*bin_spectrum_3)/sum(bin_spectrum_3).*sum(mu_3.*bin_spectrum_3)/sum(bin_spectrum_3)+N_4.*sum(mu_2.*bin_spectrum_4)/sum(bin_spectrum_4).*sum(mu_3.*bin_spectrum_4)/sum(bin_spectrum_4);
    M(3,1)=N_1.*sum(mu_3.*bin_spectrum_1)/sum(bin_spectrum_1).*sum(mu_1.*bin_spectrum_1)/sum(bin_spectrum_1)+N_2.*sum(mu_3.*bin_spectrum_2)/sum(bin_spectrum_2).*sum(mu_1.*bin_spectrum_2)/sum(bin_spectrum_2)+N_3.*sum(mu_3.*bin_spectrum_3)/sum(bin_spectrum_3).*sum(mu_1.*bin_spectrum_3)/sum(bin_spectrum_3)+N_4.*sum(mu_3.*bin_spectrum_4)/sum(bin_spectrum_4).*sum(mu_1.*bin_spectrum_4)/sum(bin_spectrum_4);
    M(3,2)=N_1.*sum(mu_3.*bin_spectrum_1)/sum(bin_spectrum_1).*sum(mu_2.*bin_spectrum_1)/sum(bin_spectrum_1)+N_2.*sum(mu_3.*bin_spectrum_2)/sum(bin_spectrum_2).*sum(mu_2.*bin_spectrum_2)/sum(bin_spectrum_2)+N_3.*sum(mu_3.*bin_spectrum_3)/sum(bin_spectrum_3).*sum(mu_2.*bin_spectrum_3)/sum(bin_spectrum_3)+N_4.*sum(mu_3.*bin_spectrum_4)/sum(bin_spectrum_4).*sum(mu_2.*bin_spectrum_4)/sum(bin_spectrum_4);
    M(3,3)=N_1.*sum(mu_3.*bin_spectrum_1)/sum(bin_spectrum_1).*sum(mu_3.*bin_spectrum_1)/sum(bin_spectrum_1)+N_2.*sum(mu_3.*bin_spectrum_2)/sum(bin_spectrum_2).*sum(mu_3.*bin_spectrum_2)/sum(bin_spectrum_2)+N_3.*sum(mu_3.*bin_spectrum_3)/sum(bin_spectrum_3).*sum(mu_3.*bin_spectrum_3)/sum(bin_spectrum_3)+N_4.*sum(mu_3.*bin_spectrum_4)/sum(bin_spectrum_4).*sum(mu_3.*bin_spectrum_4)/sum(bin_spectrum_4);
    
    cov_matrix=inv(M);
    CRLB_I(index,1) = cov_matrix(1,1);
    CRLB_Gd(index,1) = cov_matrix(2,2);
    CRLB_PMMA(index,1) = cov_matrix(3,3);
    keV(index,1) = PCD_thresh;
    
    index = index+1;
end