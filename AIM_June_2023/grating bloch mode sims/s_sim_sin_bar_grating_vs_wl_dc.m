% bloch mode simulation of sopa bar gratings
% vs. wavelength and duty cycle
% including silicon substrate

clear; close all;

% dependencies:
% addpath( genpath( 'C:\Users\bzhang\git\grating_synthesis' ) );
addpath( 'C:\Users\bz\git\grating_synthesis\main' );
addpath( 'C:\Users\bz\git\grating_synthesis\auxiliary_functions' );
addpath( 'C:\Users\beezy\git\grating_synthesis\main' );
addpath( 'C:\Users\beezy\git\grating_synthesis\auxiliary_functions' );

% save filepath
save_fpath  = 'D:\Google Drive\research\popovic group\projects\lidar SCALABLE\data\2022 04 22 - sin bar alpha dir sweep vs wl dc';
curtime     = datestr( now, 'yymmddHHMM' );

% min feat
min_width = 150;
min_space = 100;

% % pick grating
grating_type = 'lower sin bars';
lambda      = 1450:10:1650;
pitch       = 530;
duty_cycle  = fliplr( ( min_width : 10 : pitch-min_space )/pitch );
guessk      = 2.5 * 2 * pi /lambda(1);
OPTS        = struct();

% grating_type = 'upper sin bars';
% lambda      = 1450:10:1650;
% pitch       = 550;
% duty_cycle  = fliplr( ( min_width : 10 : pitch-min_space )/pitch );
% guessk      = 2.5 * 2 * pi /lambda(1);
% OPTS        = struct();

% grating_type = 'custom oxide thick sin bars';
% lambda      = 1550;
% pitch       = 530;
% duty_cycle  = 0.05:0.05:0.95;
% guessk      = 2.5 * 2 * pi /lambda;
% OPTS.custom_oxide_thickness = 300;

% grating_type = 'upper sin wg';
% lambda      = 1550;
% pitch       = 800;
% duty_cycle  = fliplr(0.05:0.05:0.95);
% % duty_cycle  = 0.05;
% % guessk      = 2.0 * 2 * pi /lambda;
% guessk      = 0.0066;

% grating_type = 'dual sin wg';
% lambda      = 1550;
% pitch       = 880;
% duty_cycle  = 0.05:0.05:0.95;
% duty_cycle  = 0.05;
% guessk      = 2.0 * 2 * pi /lambda;

% variables to save
% field_alpha = zeros( length(lambda), length(duty_cycle) ); % alpha vs. wl and duty cycle
p_rad_up    = zeros( length(lambda), length(duty_cycle) ); % pwoer radiated upwards vs. wl and duty cycle
p_rad_down  = zeros( length(lambda), length(duty_cycle) ); % power radiated downwards vs. wl and duty cycle
k           = zeros( length(lambda), length(duty_cycle) );
angle_up    = zeros( length(lambda), length(duty_cycle) );
angle_down  = zeros( length(lambda), length(duty_cycle) );

tic;
for i_wl = 1:length(lambda)
    
    fprintf('wavelength %i of %i\n', i_wl, length(lambda) );
    
    for i_dutycycle = 1:length(duty_cycle)

        fprintf('duty cycle %i of %i\n', i_dutycycle, length(duty_cycle) );

        [GC] = f_run_grating_bloch_sim(lambda(i_wl), pitch, ...
                                       grating_type, duty_cycle(i_dutycycle), guessk, OPTS );

        % update guess k and mode
        guessk                  = GC.k;
        OPTS.mode_to_overlap    = GC.E_z_for_overlap;

        % save variables
%         field_alpha(i_wl, i_dutycycle)  = imag(GC.k);
        p_rad_up(i_wl, i_dutycycle)     = GC.P_rad_up;
        p_rad_down(i_wl, i_dutycycle)   = GC.P_rad_down;
        k(i_wl, i_dutycycle)            = GC.k;
        angle_up(i_wl, i_dutycycle)     = GC.max_angle_up;
        angle_down(i_wl, i_dutycycle)   = GC.max_angle_down;
        
        % save this GC if its the first in the duty cycle loop
        if i_dutycycle == 1
            GC1 = GC;
        end
        
        toc;

    end % end duty cycle loop
    
    % update the guess k and overlap
    OPTS.mode_to_overlap    = GC1.E_z_for_overlap;
    guessk                  = GC1.k;

end

% calculate 1/e^2 power decay
decay_len_nm = 1./(2*imag(k));
decay_len_mm = 1e-6 * decay_len_nm;

% calc directionality
pdir = p_rad_up./p_rad_down;

% save the data
save_fname = [ curtime '_pitch' num2str(pitch) '_' strrep(grating_type, ' ', '_') ];
save( [ save_fpath filesep save_fname '.mat' ] );

% plot 1/e^2 power decay length as function of wl and duty cycle
figure;
imagesc( duty_cycle, lambda, decay_len_mm ); hold on;
contour( duty_cycle, lambda, decay_len_mm, 'w', 'showtext', 'on' );
xlabel('duty cycle'); ylabel('wavelength (nm)');
title('1/e^2 power decay length');
colorbar;
colormap('hsv');
set(gca, 'ydir', 'normal');

% plot power directionality 
figure;
imagesc( duty_cycle, lambda, pdir ); hold on;
contour( duty_cycle, lambda, pdir, 'w', 'showtext', 'on' );
xlabel('duty cycle'); ylabel('wavelength (nm)');
title('power unidirectionality (P up/P down)');
colorbar;
set(gca, 'ydir', 'normal');

% plot radiating angle
figure;
imagesc( duty_cycle, lambda, angle_up ); hold on;
contour( duty_cycle, lambda, angle_up, 'w', 'showtext', 'on' );
xlabel('duty cycle'); ylabel('wavelength (nm)');
title('radiating angle');
colorbar;
set(gca, 'ydir', 'normal');

% calculate insertion loss for given length of grating
n_rows          = 96;
l_gcrow         = 1.8; % mm
l_total         = n_rows * l_gcrow; % mm
prad            = 1 - exp( - 2 * imag(k) * 1e6 * l_total ); % total radiated power
pradloss_dB     = -10*log10(prad);
pdirloss_dB     = -10*log10(1./(1 + 1./pdir)); % directionality loss
rad_up_IL_dB    = pradloss_dB + pdirloss_dB;

% plot insertion loss contributed by power directionality
figure;
imagesc( duty_cycle, lambda, pdirloss_dB ); hold on;
contour( duty_cycle, lambda, pdirloss_dB, 'w', 'showtext', 'on' );
xlabel('duty cycle'); ylabel('wavelength (nm)');
title('Insertion loss (dB) contrib. by power directionality');
colorbar;
colormap('jet');
set(gca, 'ydir', 'normal');

% plot insertion loss contributed by total radiated power
figure;
imagesc( duty_cycle, lambda, pradloss_dB ); hold on;
contour( duty_cycle, lambda, pradloss_dB, 'w', 'showtext', 'on' );
xlabel('duty cycle'); ylabel('wavelength (nm)');
title('Insertion loss (dB) contrib. by total radiated power');
colorbar;
colormap('jet');
set(gca, 'ydir', 'normal');

% plot insertion loss
figure;
imagesc( duty_cycle, lambda, rad_up_IL_dB ); hold on;
contour( duty_cycle, lambda, rad_up_IL_dB, 'w', 'showtext', 'on' );
xlabel('duty cycle'); ylabel('wavelength (nm)');
title('Upwards radiation insertion loss (dB)');
colorbar;
colormap('jet');
text( 0.7, 1425, ['n rows ' num2str(n_rows) ', l row ' num2str(l_gcrow) ' mm'], 'fontsize', 12);
set(gca, 'ydir', 'normal');


% plot these things at 50% duty cycle vs lambda
interp2_dc50        = @( x ) interp2( duty_cycle, lambda, x, 0.5, lambda );
decay_len_mm_dc50   = interp2_dc50( decay_len_mm );
pdir_dc50           = interp2_dc50( pdir );
angle_up_dc50       = interp2_dc50( angle_up );
pradloss_dB_dc50    = interp2_dc50( pradloss_dB );
pdirloss_dB_dc50    = interp2_dc50( pdirloss_dB );
rad_up_IL_dB_dc50   = interp2_dc50( rad_up_IL_dB );

figure;
plot( lambda, decay_len_mm_dc50, '-o' );
xlabel('wavelength'); ylabel('1/e^2 decay length, mm');
makeFigureNice();

figure;
plot( lambda, pdir_dc50, '-o' );
xlabel('wavelength'); ylabel('P up/P down');
makeFigureNice();

figure;
plot( lambda, angle_up_dc50, '-o' );
xlabel('wavelength'); ylabel('upwards radiation angle');
makeFigureNice();

figure;
plot( lambda, pradloss_dB_dc50, '-o' ); hold on;
plot( lambda, pdirloss_dB_dc50, '-o' );
plot( lambda, rad_up_IL_dB_dc50, '-o' );
xlabel('wavelength'); ylabel('insertion loss (dB)');
title(['n rows ' num2str(n_rows) ', l row ' num2str(l_gcrow) ' mm, 50% dutycycle']);
makeFigureNice();


% % calculelate 1/e distance
% decay_n_nm = 1./(2*field_alpha);
% decay_len_mm = 1e-6 * decay_len_nm;
% decay_len_cm = 1e-7 * decay_len_nm;
%     
% % calc loss in dB/cm
% loss_dbcm = -10*log10( exp( -2 * field_alpha * 1e7 ) );
% loss_dbmm = loss_dbcm/10;
% 
% % calculate loss for given sopa length
% sopa_grating_length_cm  = 2;
% loss_total              = loss_dbcm * sopa_grating_length_cm;
% 
% % calculate radiation for given sopa length
% radiated_power_db = 10*log10( 1 - exp( -2 * field_alpha * 1e7 * sopa_grating_length_cm ) );

% % plot decay length vs. duty cycle
% figure;
% plot( duty_cycle, decay_len_cm, '-o' );
% xlabel('duty cycle'); ylabel('1/e power decay length (cm)');
% title( grating_type );
% makeFigureNice();

% % plot loss vs. duty cycle
% figure;
% plot( duty_cycle, loss_dbmm, '-o' );
% xlabel('duty cycle'); ylabel('radiated power in dB/mm');
% title( grating_type );
% makeFigureNice();
% 
% % plot loss vs. grating longitudinal width
% figure;
% plot( duty_cycle.*pitch, loss_dbmm, '-o' );
% xlabel('grating longitudinal width (nm)'); ylabel('radiated power in dB/mm');
% title( grating_type );
% makeFigureNice();

% % plot loss vs. duty cycle
% figure;
% plot( duty_cycle, loss_total, '-o' );
% xlabel('duty cycle'); ylabel('loss in dB');
% title( [ grating_type ' ' num2str(sopa_grating_length_cm) ' cm grating' ] );
% makeFigureNice();

% % plot radiated power vs. duty cycle
% figure;
% plot( duty_cycle, radiated_power_db, '-o' );
% xlabel('duty cycle'); ylabel('radiated power in dB');
% title( [ grating_type ' ' num2str(sopa_grating_length_cm) ' cm grating' ] );
% makeFigureNice();

% 
% % plot the field
% figure;
% imagesc( GC.x_coords, GC.y_coords, real(GC.E_z) );
% xlabel('x');
% ylabel('y');
% colorbar;
% colormap('redbluehilight');
% set(gca,'ydir','normal');
% caxis( [-1e-4, 1e-4 ] );











