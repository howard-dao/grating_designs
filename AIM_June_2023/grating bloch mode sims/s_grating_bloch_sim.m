% bloch mode simulation of sopa bar gratings
% vs. type and duty cycle

clear; close all;

% dependencies:
% addpath( genpath( 'C:\Users\bzhang\git\grating_synthesis' ) );
addpath( 'C:\Users\bz\git\grating_synthesis\main' );
addpath( 'C:\Users\bz\git\grating_synthesis\auxiliary_functions' );
addpath( 'C:\Users\beezy\git\grating_synthesis\main' );
addpath( 'C:\Users\beezy\git\grating_synthesis\auxiliary_functions' );


% % pick grating
% grating_type = 'lower sin bars';
% lambda      = 1550;
% pitch       = 530;
% duty_cycle  = 0.05:0.05:0.95;
% guessk      = 2.5 * 2 * pi /lambda;

grating_type = 'upper sin bars';
lambda      = 1550;
pitch       = 530;
duty_cycle  = 0.8; %0.05:0.05:0.95;
guessk      = 2.5 * 2 * pi /lambda;
OPTS        = struct();

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
field_alpha = zeros( size(duty_cycle) );

for i_dutycycle = 1:length(duty_cycle)

    fprintf('duty cycle %i of %i\n', i_dutycycle, length(duty_cycle) );
    
    if i_dutycycle > 1
        OPTS.mode_to_overlap = GC.E_z_for_overlap;
    end
    
    tic;
    [GC] = f_run_grating_bloch_sim(lambda, pitch, ...
                                           grating_type, duty_cycle(i_dutycycle), guessk, OPTS );
    toc;
    
    % update guess k
    guessk = GC.k;
    
    % save alpha
    field_alpha(i_dutycycle) = imag(GC.k);

end

    
% calculate 1/e distance
decay_len_nm = 1./(2*field_alpha);
decay_len_mm = 1e-6 * decay_len_nm;
decay_len_cm = 1e-7 * decay_len_nm;
    
% calc loss in dB/cm
loss_dbcm = -10*log10( exp( -2 * field_alpha * 1e7 ) );
loss_dbmm = loss_dbcm/10;

% calculate loss for given sopa length
sopa_grating_length_cm  = 2;
loss_total              = loss_dbcm * sopa_grating_length_cm;

% calculate radiation for given sopa length
radiated_power_db = 10*log10( 1 - exp( -2 * field_alpha * 1e7 * sopa_grating_length_cm ) );

% % plot decay length vs. duty cycle
% figure;
% plot( duty_cycle, decay_len_cm, '-o' );
% xlabel('duty cycle'); ylabel('1/e power decay length (cm)');
% title( grating_type );
% makeFigureNice();

% plot loss vs. duty cycle
figure;
plot( duty_cycle, loss_dbmm, '-o' );
xlabel('duty cycle'); ylabel('radiated power in dB/mm');
title( grating_type );
makeFigureNice();

% plot loss vs. grating longitudinal width
figure;
plot( duty_cycle.*pitch, loss_dbmm, '-o' );
xlabel('grating longitudinal width (nm)'); ylabel('radiated power in dB/mm');
title( grating_type );
makeFigureNice();

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











