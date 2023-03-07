% loading and plotting things im interested in for bar gratings

clear; close all;

fig_pos = [300, 230,400, 286];

% data
filepath = 'G:\My Drive\research\popovic group\projects\lidar SCALABLE\data\2022 04 22 - sin bar alpha dir sweep vs wl dc';
filename_uppersin = '2204221322_pitch530_upper_sin_bars.mat';
filename_lowersin = '2204221556_pitch530_lower_sin_bars.mat';

% grating size
n_rows          = 48;
l_gcrow         = 0.8; % mm
l_total         = n_rows * l_gcrow; % mm

% load and plot lower sin
lowersin_data = load([ filepath filesep filename_lowersin]);

% plot 1/e^2 power decay length as function of wl and duty cycle
figure('name', 'decaylength-lowersin', 'position', fig_pos);
imagesc( lowersin_data.duty_cycle, lowersin_data.lambda, lowersin_data.decay_len_mm ); hold on;
contour( lowersin_data.duty_cycle, lowersin_data.lambda, lowersin_data.decay_len_mm, 'k', 'showtext', 'on' );
xlabel('Duty cycle'); ylabel('Wavelength (nm)');
title('1/e^2 power decay length (mm)');
colorbar;
colormap(hsv);
set(gca, 'ydir', 'normal');

% plot power directionality 
figure('name', 'directionality-lowersin', 'position', fig_pos);
imagesc( lowersin_data.duty_cycle, lowersin_data.lambda, lowersin_data.pdir ); hold on;
% contour( lowersin_data.duty_cycle, lowersin_data.lambda, lowersin_data.pdir, 0:4, 'w', 'showtext', 'on' );
[C,h] = contour( lowersin_data.duty_cycle, lowersin_data.lambda, lowersin_data.pdir, 0:4, 'k', 'showtext', 'on' );
% clabel(C,h,'FontSize',15,'Color','r','FontName','Courier')
clabel(C,h,'FontSize',10,'Color','k');
xlabel('Duty cycle'); ylabel('Wavelength (nm)');
title('Power unidirectionality (P_{up}/P_{down})');
colorbar;
colormap(jet);
set(gca, 'ydir', 'normal');

% calculate insertion loss for given length of grating
prad            = 1 - exp( - 2 * imag(lowersin_data.k) * 1e6 * l_total ); % total radiated power
pradloss_dB     = -10*log10(prad);
pdirloss_dB     = -10*log10(1./(1 + 1./lowersin_data.pdir)); % directionality loss
rad_up_IL_dB    = pradloss_dB + pdirloss_dB;

% plot insertion loss
figure('name', 'IL-lowersin', 'position', fig_pos);
imagesc( lowersin_data.duty_cycle, lowersin_data.lambda, rad_up_IL_dB ); hold on;
[C,h] = contour( lowersin_data.duty_cycle, lowersin_data.lambda, rad_up_IL_dB, 1:4, 'k', 'showtext', 'on' );
clabel(C,h,'FontSize',10,'Color','k');
xlabel('Duty cycle'); ylabel('Wavelength (nm)');
title('Insertion loss (dB)');
colorbar;
colormap(jet);
% text( 0.7, 1425, ['n rows ' num2str(n_rows) ', l row ' num2str(l_gcrow) ' mm'], 'fontsize', 12);
set(gca, 'ydir', 'normal');


% load and plot lower sin
uppersin_data = load([ filepath filesep filename_uppersin]);

% plot 1/e^2 power decay length as function of wl and duty cycle
figure('name', 'decaylength-uppersin', 'position', fig_pos);
imagesc( uppersin_data.duty_cycle, uppersin_data.lambda, 0.1*uppersin_data.decay_len_mm ); hold on;
[C,h] = contour( uppersin_data.duty_cycle, uppersin_data.lambda, 0.1*uppersin_data.decay_len_mm, [25,50,100,200,300], 'k', 'showtext', 'on' );
clabel(C,h,'FontSize',10,'Color','k');
xlabel('Duty cycle'); ylabel('Wavelength (nm)');
title('1/e^2 power decay length (cm)');
colorbar;
colormap(hsv);
set(gca, 'ydir', 'normal');

% plot power directionality 
figure('name', 'directionality-uppersin', 'position', fig_pos);
imagesc( uppersin_data.duty_cycle, uppersin_data.lambda, uppersin_data.pdir ); hold on;
% contour( lowersin_data.duty_cycle, lowersin_data.lambda, lowersin_data.pdir, 0:4, 'w', 'showtext', 'on' );
[C,h] = contour( uppersin_data.duty_cycle, uppersin_data.lambda, uppersin_data.pdir, 0:4, 'k', 'showtext', 'on' );
% clabel(C,h,'FontSize',15,'Color','r','FontName','Courier')
clabel(C,h,'FontSize',10,'Color','k');
xlabel('Duty cycle'); ylabel('Wavelength (nm)');
title('Power unidirectionality (P_{up}/P_{down})');
colorbar;
colormap(jet);
set(gca, 'ydir', 'normal');

% calculate insertion loss for given length of grating
prad            = 1 - exp( - 2 * imag(uppersin_data.k) * 1e6 * l_total ); % total radiated power
pradloss_dB     = -10*log10(prad);
pdirloss_dB     = -10*log10(1./(1 + 1./uppersin_data.pdir)); % directionality loss
rad_up_IL_dB    = pradloss_dB + pdirloss_dB;

% plot insertion loss
figure('name', 'IL-uppersin', 'position', fig_pos);
imagesc( uppersin_data.duty_cycle, uppersin_data.lambda, rad_up_IL_dB ); hold on;
[C,h] = contour( uppersin_data.duty_cycle, uppersin_data.lambda, rad_up_IL_dB, 10:2:20, 'k', 'showtext', 'on' );
clabel(C,h,'FontSize',10,'Color','k');
xlabel('Duty cycle'); ylabel('Wavelength (nm)');
title('Insertion loss (dB)');
colorbar;
colormap(jet);
% text( 0.7, 1425, ['n rows ' num2str(n_rows) ', l row ' num2str(l_gcrow) ' mm'], 'fontsize', 12);
set(gca, 'ydir', 'normal');

% save figures
save_fig_path = 'G:\My Drive\research\popovic group\thesis\figures\by section\7 3 2\working';
save_all_figs(save_fig_path);






















