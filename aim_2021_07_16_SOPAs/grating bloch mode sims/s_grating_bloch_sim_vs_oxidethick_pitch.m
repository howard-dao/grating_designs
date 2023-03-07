% bloch mode simulation of sopa bar gratings
% vs. custom oxide thickness and duty cycle

clear; close all;

% dependencies:
% addpath( genpath( 'C:\Users\bzhang\git\grating_synthesis' ) );
addpath( 'C:\Users\bz\git\grating_synthesis\main' );
addpath( 'C:\Users\bz\git\grating_synthesis\auxiliary_functions' );
addpath( 'C:\Users\beezy\git\grating_synthesis\main' );
addpath( 'C:\Users\beezy\git\grating_synthesis\auxiliary_functions' );

% grating settings
grating_type = 'custom oxide thick sin bars';
lambda      = 1550;
pitch       = 530;
duty_cycle  = 0.05:0.05:0.95;
guessk      = 2.5 * 2 * pi /lambda;
oxide_thicknesses = 100 : 50 : 420;

% variables to save
field_alpha     = zeros( length(oxide_thicknesses), length(duty_cycle) ); % dimensions oxide thick vs. duty cycle
rad_angle       = zeros( length(oxide_thicknesses), length(duty_cycle) ); % dimensions oxide thick vs. duty cycle
directionality  = zeros( length(oxide_thicknesses), length(duty_cycle) ); % power up/power down dimensions oxide thick vs. duty cycle

iloop = 1;
for i_oxthick = 1:length(oxide_thicknesses)
    
    fprintf('loop %i of %i\n', iloop, length(oxide_thicknesses)*length(duty_cycle));
    OPTS.custom_oxide_thickness = oxide_thicknesses(i_oxthick);
    
    for i_dutycycle = 1:length(duty_cycle)

        if i_dutycycle > 1
            OPTS.mode_to_overlap = GC.E_z_for_overlap;
        end

        tic;
        [GC] = f_run_grating_bloch_sim(lambda, pitch, ...
                                               grating_type, duty_cycle(i_dutycycle), guessk, OPTS );
        toc;

        % update guess k
        guessk = GC.k;

        % save variables
        field_alpha(i_oxthick,i_dutycycle)      = imag(GC.k);
        rad_angle(i_oxthick,i_dutycycle)        = GC.max_angle_up;
        directionality(i_oxthick,i_dutycycle)   = GC.directivity;
        
        iloop  = iloop + 1;

    end
    
    guessk = 2.5 * 2 * pi /lambda;
    
end
    
% calculate 1/e distance
decay_len_nm = 1./(2*field_alpha);
decay_len_mm = 1e-6 * decay_len_nm;
decay_len_cm = 1e-7 * decay_len_nm;
    
% calc loss in dB/cm
loss_dbcm = -10*log10( exp( -2 * field_alpha * 1e7 ) );
loss_dbmm = loss_dbcm/10;

% plot loss vs. oxide thickness and duty cycle
figure;
imagesc( duty_cycle, oxide_thicknesses, loss_dbmm );
xlabel('duty cycle'); ylabel('oxide thickness (nm)');
title( 'Radiation loss (dB/mm)' );
colorbar;
set(gca, 'ydir', 'normal');

% plot loss vs. oxide thickness and duty cycle
figure;
imagesc( duty_cycle.*pitch, oxide_thicknesses, loss_dbmm );
xlabel('grating longitudinal width (nm)'); ylabel('oxide thickness (nm)');
title( 'Radiation loss (dB/mm)' );
colorbar;
set(gca, 'ydir', 'normal');

% plot radiation angle vs. oxide thickness and duty cycle
figure;
imagesc( duty_cycle.*pitch, oxide_thicknesses, rad_angle );
xlabel('grating longitudinal width (nm)'); ylabel('oxide thickness (nm)');
title( 'Radiation angle (deg)' );
colorbar;
set(gca, 'ydir', 'normal');

% plot directionality vs. oxide thickness and duty cycle
figure;
imagesc( duty_cycle.*pitch, oxide_thicknesses, directionality );
xlabel('grating longitudinal width (nm)'); ylabel('oxide thickness (nm)');
title( 'Directionality (up/down)' );
colorbar;
set(gca, 'ydir', 'normal');

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











