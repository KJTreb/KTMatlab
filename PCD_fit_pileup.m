%% Generate LUT: Input = [PCD counts, EID signal]; Output = N_PCD_corr
close all; clear; clc;


mA = [2, 10, 16, 20, 25, 32, 40, 50, 63, 80, 100, 125, 160, 200];
N = numel(mA);
x=mA';

rows = 64;
cols = 5120;

prj_TE = zeros(rows, cols, N);

for i = 1:N
    filename = sprintf('data/%dmA_HE.raw', mA(i));
    temp = MgReadRawFile(filename, rows, cols, 50, 0, 0, 'float32');
    prj_TE(:,:,i) = mean(temp,3);
    clear temp
end

p_PCD_HE = cell(rows, cols);
%argmax_mA = cell(rows, cols);
Rsquare = cell(rows, cols);

%dead time model fit
%fitfun = fittype('w.*(a.*x.*exp(-b.*a.*x))+(1-w).*(a.*x/(1+b.*a.*x))', 'dependent', {'y'}, 'independent', {'x'}, 'coefficients', {'a', 'b', 'w'} );
%options = fitoptions(fitfun);
% options.Lower = [0 0 0];
% options.Upper = [Inf Inf 1];
% options.StartPoint = [1 1 0.5];
fitfun = fittype('p1.*x.^5 + p2.*x.^4 + p3.*x.^3 + p4.*x.^2 + p5.*x + p6', 'dependent', {'y'}, 'independent', {'x'}, 'coefficients', {'p1', 'p2', 'p3', 'p4', 'p5', 'p6'} );
options = fitoptions(fitfun);
options.StartPoint = [1 1 1 1 22.81 0];
options.Lower = [-Inf -Inf -Inf -Inf 22.55 0];
options.Upper = [Inf Inf Inf Inf 23.07, 0];

pb = MgCmdLineProgressBar('Calculating: ');
tic
for col = 1:cols
    %pb.print(col, cols);
    for row = 1:rows
        % get projection data
        % TE
        p_TE = prj_TE(row, col, :);
        p_TE = squeeze(p_TE)';

        y = p_TE';
        %fit to model
        [fitted_curve,gof] = fit(x,y,fitfun, options);

        % get fit parameters
        p_PCD_HE{row, col} = coeffvalues(fitted_curve);
        Rsquare{row, col} = gof.rsquare;
    end
    time_elapsed = toc/60;
    time_remaining = (toc/col)*(cols-col)/60;
    clc
    disp(['Time elapsed: ',num2str(time_elapsed),' minutes. Time remaining: ',num2str(time_remaining),' minutes.'])
    disp([num2str(col), '/', num2str(cols)])
end

MgMkdir('coefficients', false);
save('coefficients/p_PCD_HE_polynomial', 'p_PCD_HE');
%save('coefficients/argmax_mA', 'argmax_mA')