% bloch mode simulation of sopa bar gratings
% vs. wavelength

clear; close all;

% dependencies:
% addpath( genpath( 'C:\Users\bzhang\git\grating_synthesis' ) );
addpath( 'C:\Users\bz\git\grating_synthesis\main' );
addpath( 'C:\Users\bz\git\grating_synthesis\auxiliary_functions' );
addpath( 'C:\Users\beezy\git\grating_synthesis\main' );
addpath( 'C:\Users\beezy\git\grating_synthesis\auxiliary_functions' );

% grating settings
grating_type = 'custom oxide thick sin bars';
lambda0     = 1550;
pitch       = 530;
% duty_cycle  = 0.05:0.05:0.95;
duty_cycle  = 0.5;
guess_neff  = 2.9;
% guessk      = 2.5 * 2 * pi /lambda0;
oxide_thickness = 200;

% wavelengths to sweep
lambdas = fliplr( (lambda0 - 100) : 50 : (lambda0 + 100) );
% lambdas = 1550 : 50 : 1650;

% variables to save
field_alpha     = zeros( size(lambdas) ); % dimensions oxide thick vs. duty cycle
rad_angle       = zeros( size(lambdas) ); % dimensions oxide thick vs. duty cycle
directionality  = zeros( size(lambdas) ); % power up/power down dimensions oxide thick vs. duty cycle

% options
OPTS.custom_oxide_thickness = oxide_thickness;

guessk = guess_neff * 2*pi/lambdas(1);

iloop = 1;
for i_lambda = 1:length(lambdas)
    
    fprintf('loop %i of %i\n', i_lambda, length(lambdas));
    
    if i_lambda > 1
        OPTS.mode_to_overlap = GC.E_z_for_overlap;
    end

    tic;
    [GC] = f_run_grating_bloch_sim(lambdas(i_lambda), pitch, ...
                                           grating_type, duty_cycle, guessk, OPTS );
    toc;

    % update guess k
    guessk = GC.k;

    % save variables
    field_alpha(i_lambda)      = imag(GC.k);
    rad_angle(i_lambda)        = GC.max_angle_up;
    directionality(i_lambda)   = GC.directivity;
    
end
    
% calculate 1/e distance
decay_len_nm = 1./(2*field_alpha);
decay_len_mm = 1e-6 * decay_len_nm;
decay_len_cm = 1e-7 * decay_len_nm;
    
% calc loss in dB/cm
loss_dbcm = -10*log10( exp( -2 * field_alpha * 1e7 ) );
loss_dbmm = loss_dbcm/10;

% plot loss vs. wavelength
figure;
plot(lambdas, loss_dbmm, '-o' );
xlabel('wavelength'); ylabel('radiation loss (db/mm)');
title(['period ' num2str(pitch) ', duty cycle ' num2str(duty_cycle) ', oxide thick ' num2str(oxide_thickness) ]);
makeFigureNice();

% plot angle vs. wavelength
figure;
plot(lambdas, rad_angle, '-o' );
xlabel('wavelength'); ylabel('radiation angle (degrees)');
title(['period ' num2str(pitch) ', duty cycle ' num2str(duty_cycle) ', oxide thick ' num2str(oxide_thickness) ]);
makeFigureNice();

% plot directionality vs. wavelength
figure;
plot(lambdas, directionality, '-o' );
xlabel('wavelength'); ylabel('directionality (up/down)');
title(['period ' num2str(pitch) ', duty cycle ' num2str(duty_cycle) ', oxide thick ' num2str(oxide_thickness) ]);
makeFigureNice();

















