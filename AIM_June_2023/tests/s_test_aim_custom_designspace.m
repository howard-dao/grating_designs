% custom aim gc cell simulation for SOPAIPILLA
% test the functionality of design space code

clear; close all;

% dependencies
addpath(genpath('C:\Users\howar\Siphot\Grating Couplers\grating_synthesis'));
addpath(genpath('C:\Users\howar\Siphot\Grating Couplers\grating_designs\AIM_June_2023'));

disc = 10;
lambda = 1550;
optimal_angle = 0;
coupling_direction = 'up';

% initial settings
units               = 'nm';
background_index    = 1.445;
y_domain_size       = 4e3;
data_notes          = [ 'lambda ', num2str(lambda), ...
                        ' optimal angle ', num2str(optimal_angle), ...
                        ' discretization ', num2str(disc), ...
                        ' coupling direction ', coupling_direction];

% display inputs
fprintf('Running design space synthesis on aim grating\n');
fprintf('Inputs are:\n');
fprintf('Wavelength: %f %s\n', lambda, units);
fprintf('Angle: %f degrees\n', optimal_angle);
fprintf('Coupling direction: %s\n', coupling_direction);
fprintf('Discretization: %f\n', disc);
fprintf('Designed for coupling into: %s\n', 'oxide index 1.445');

% geometry parameters
thick_Si = 220;
space_Si_SiN = 100;
thick_bot_SiN = 60;
space_SiN_SiN = 0;
thick_top_SiN = 160;

geometry = struct('thick_Si', thick_Si, ...
                    'space_Si_SiN', space_Si_SiN, ...
                    'thick_bot_SiN', thick_bot_SiN, ...
                    'space_SiN_SiN', space_SiN_SiN, ...
                    'thick_top_SiN', thick_top_SiN ...
                    );

h_makeGratingCell = @(dxy, background_index, y_domain_size, period, fill_top, fill_bot, offset_ratio) ...
                        f_makeGratingCell_AIM_custom(   dxy, ...
                                                        lambda, ...
                                                        y_domain_size, ...
                                                        period, ...
                                                        fill_top, ...
                                                        fill_bot, ...
                                                        offset_ratio, ...
                                                        geometry ...
                                                        );

% make synthesis object
synth_obj = c_synthTwoLevelGrating( 'discretization',    disc, ...
                                    'units',             units, ...
                                    'lambda',            lambda, ...
                                    'background_index',  background_index, ...
                                    'y_domain_size',     y_domain_size, ...
                                    'optimal_angle',     optimal_angle, ...
                                    'data_notes',        data_notes, ...
                                    'coupling_direction', coupling_direction, ...
                                    'h_makeGratingCell', h_makeGratingCell ...
                                    );


% % grating cell settings
% dxy             = 10;
% lambda          = 1550;
% k0              = 2*pi/lambda;
% y_domain_size   = 4e3;
% 
% neff_slab   = 2.842431701;
% n_clad      = 1.45;
% theta       = 0; % deg
% period      = 2*pi/( k0 * ( neff_slab - n_clad*sin( theta*pi/180 ) ) );
% period      = round(period/dxy)  * dxy;

% duty cycles
dc_bots = 0.48 : 0.02 : 0.50;
dc_tops = 0.48 : 0.02 : 0.50;

% % run waveguide simulation for initial guesses
% waveguide = f_makeGratingCell_AIM_custom(dxy, lambda, y_domain_size, 2*dxy, 1.0, 1.0, 0.0, geometry);
% 
% guess_n = 1.0 * max(waveguide.N(:));
% guess_k = guess_n * 2*pi/lambda;
% 
% num_wg_modes = 5;
% BC = 0;
% pml_options_wg = [0, 200, 20, 2];
% 
% waveguide = waveguide.runSimulation(num_wg_modes, BC, pml_options_wg, k0, guess_k);
% 
% guess_k             = waveguide.k;
% guess_offset        = 0;
% n_optimize_loops    = 2;
% 
% % calculate analytical period which would approximately phase-match to
% % desired output angle
% waveguide_k     = waveguide.k;
% k0              = waveguide.background_index * 2*pi/lambda;
% kx              = k0 * sin( (pi/180) * theta);
% guess_period    = 2*pi / (waveguide_k - kx);
% guess_period    = round(guess_period / dxy) * dxy;
% 
% [waveguide, e_z_overlap_ext] = waveguide.stitch_E_field(waveguide.Phi, ...
%                                                         real(waveguide.k), ...
%                                                         round(guess_period/waveguide.domain_size(2)) ...
%                                                         );
% waveguide.E_z_for_overlap = e_z_overlap_ext;
% 
% % set grating solver settings
% num_modes   = 5;
% BC          = 0;    % 0 = PEC
% pml_options = [1, 100, 20, 2];
% OPTS        = struct( 'mode_to_overlap', e_z_overlap_ext );
% sim_opts    = struct('num_modes', num_modes, 'BC', BC, 'pml_options', pml_options);
% 
% guess_GC = waveguide;
% 
% % initialize variables
% directivities_vs_fills  = zeros(length(dc_bots), length(dc_tops));
% angles_vs_fills         = zeros(length(dc_bots), length(dc_tops));
% periods_vs_fills        = zeros(length(dc_bots), length(dc_tops));
% offsets_vs_fills        = zeros(length(dc_bots), length(dc_tops));
% k_vs_fills              = zeros(length(dc_bots), length(dc_tops));
% GC_vs_fills             = cell( length(dc_bots), length(dc_tops));
% scatter_str_vs_fills    = zeros(length(dc_bots), length(dc_tops));

% run design space generation, top bot version
tic;
synth_obj = synth_obj.generate_design_space(dc_bots, dc_tops);
toc;

fprintf('Design space generation sweep is done\n');
fprintf('Saving data...\n');

% make folder to save to
save_data_path = [  pwd filesep 'AIM_TEST'... 
                    '__custom_thicknesses_' num2str(thick_Si) ...
                    '_' num2str(space_Si_SiN) ...
                    '_' num2str(thick_bot_SiN) ...
                    '_' num2str(space_SiN_SiN) ...
                    '_' num2str(thick_top_SiN) ...
                    '_nm__lam' num2str(lambda) ...
                    '__ang' num2str(optimal_angle) ...
                    '__dx' num2str(synth_obj.discretization) ...
                    '_' coupling_direction ];
mkdir( save_data_path );

% clear the GC from the data and save
save( [ save_data_path filesep 'synth_obj' ], 'synth_obj' );

fprintf('Data saved\n');

% make folder to save plots to
save_plots_path = [ save_data_path filesep 'plots' ];
mkdir( save_plots_path );

fig_suffix = '';

f_save_design_space_plots( synth_obj, fig_suffix, save_plots_path );









