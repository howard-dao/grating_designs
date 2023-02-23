function [] = f_run_designspace_duallevel_45CLO( lambda, optimal_angle, shallow_or_full, coupling_direction, disc )
% authors: bohan
% 
% runs dual level unit cell generator
%
% Inputs
%   lambda
%       wavelength (nm)
%   optimal_angle
%       desired output angle
%   shallow_or_full
%       'shallow' for shallow etch, 'full' for full etch
%   coupling_direction
%       'up' or 'down'
%   disc
%       discretization (nm)

% dependencies
% all synth code
% laptop
addpath( genpath( [ 'C:\Users\howar\OneDrive\Desktop\Documents' ...
                    '\Boston University\Silicon Photonics\Grating Couplers\grating_synthesis'] ) );
% SCC
addpath( genpath( '\projectnb\siphot\howard\git\grating_synthesis' ) );

% initial settings
units               = 'nm';
index_clad          = 1.45;              % 1.448;
y_domain_size       = 4000;
data_notes          = ['lambda ' num2str(lambda) ' optimal angle ' num2str(optimal_angle) ' etch depth ' shallow_or_full ...
                       ' coupling direction ' coupling_direction ' discretization ' num2str(disc) ];

% display inputs
fprintf('Inputs are:\n');
fprintf('Wavelength: %f %s\n', lambda, units);
fprintf('Angle: %f degrees\n', optimal_angle);
fprintf('Etch depth: %s\n\n', shallow_or_full);
fprintf('Coupling direction: %s\n\n', coupling_direction);
fprintf('Discretization: %f\n\n', disc);

% handle to grating cell gen function for 45SPCLO
h_makeGratingCell   = @( dxy, background_index, y_domain_size, ...
                         period, fill_top, fill_bot, offset_ratio ) ...
                        f_makeGratingCell_45CLO( dxy, lambda, y_domain_size, ...
                                         period, fill_top, fill_bot, offset_ratio, shallow_or_full );

% make synthesis object
synth_obj = c_synthTwoLevelGrating(   'discretization',    disc, ...
                                      'units',             units,   ...
                                      'lambda',            lambda, ...
                                      'background_index',  index_clad,    ...
                                      'y_domain_size',     y_domain_size, ...
                                      'optimal_angle',     optimal_angle, ...
                                      'data_notes',        data_notes, ...
                                      'coupling_direction', coupling_direction, ...
                                      'h_makeGratingCell', h_makeGratingCell ...
                                      );

% display Q for logging purposes
synth_obj

% design space fills for top bot version
fill_bots = 0.02:0.02:0.98;
fill_tops = 0.02:0.02:0.98;

% for debugging
%fill_bots = 0.96;
%fill_tops = 0.96;

% fill_bots = 0.93:0.02:0.95;
% fill_tops = 0.93:0.02:0.95;

% fill_bots = 0.8;
% fill_tops = 0.8;


% run design space generation, top bot version
tic;
synth_obj = synth_obj.generate_design_space( fill_bots, fill_tops );
toc;

fprintf('Design space generation sweep is done\n');
fprintf('Saving data...\n');

% make folder to save to
disp(pwd)
disp(synth_obj.start_time)
disp(optimal_angle)
disp(shallow_or_full)
disp(synth_obj.discretization)
save_data_path = [ pwd filesep synth_obj.start_time ...
                    '45CLO_lam' num2str(lambda) ...
                    '_ang' num2str(optimal_angle) ...
                    '_etch_' shallow_or_full ...
                    '_dx_' num2str(synth_obj.discretization) ... 
                    '_' coupling_direction ];
mkdir( save_data_path );

% clear the GC from the data and save
% save( [ save_data_path filesep 'synth_obj_' synth_obj.start_time 'lambda' num2str(lambda) '_optangle' num2str(optimal_angle) '_etch_' shallow_or_full '_NO_GC' '.mat' ], 'synth_obj', '-v7.3' );
save( [ save_data_path filesep 'synth_obj_45CLO' ] );

fprintf('Data saved\n');

% make folder to save plots to
save_plots_path = [ save_data_path filesep 'plots' ];
mkdir( save_plots_path );

fig_suffix = ''; %[ 'lambda' num2str(lambda) '_optangle' num2str(optimal_angle) '_etch_' shallow_or_full ];

f_save_design_space_plots( synth_obj, fig_suffix, save_plots_path );


end


        
