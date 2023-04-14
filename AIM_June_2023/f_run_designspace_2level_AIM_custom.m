function [] = f_run_designspace_2level_AIM_custom(  lambda, optimal_angle, disc, coupling_direction, ...
                                                    thick_Si, space_Si_SiN, thick_bot_SiN, ...
                                                    space_SiN_SiN, thick_top_SiN)
% authors: bohan, howard
% 
% script for running the bi-level design space generation for aim layers
%
% saves results to a synth_obj
%
% Inputs
%   lambda
%       type: double, scalar
%       desc: wavelength [nm]
%   optimal_angle
%       type: double, scalar
%       desc: desired output angle [degrees]
%   disc
%       type: double, scalar
%       desc: discretization [nm]
%   coupling_direction
%       type: string
%       desc: either 'up' or 'down'
%   thick_Si
%	    type: double, scalar
%	    desc: thickness of the silicon waveguide [nm]
%   space_Si_SiN
%	    type: double, scalar
%	    desc: spacing between the silicon and bottom nitride layers [nm]
%   thick_bot_SiN
%	    type: double, scalar
%	    desc: thickness of bottom nitride [nm]
%   space_SiN_SiN
%	    type: double, scalar
%	    desc: spacing between the bottom and top nitride layers [nm]
%   thick_top_sin
%       type: double, scalar
%       desc: thickness of top nitride layer [nm]

% dependencies
% desktop
addpath(genpath('C:\Users\hdao\git\grating_synthesis'));
% laptop
addpath(genpath('C:\Users\howar\Siphot\Grating Couplers\grating_synthesis'));
% SCC
addpath(genpath('/projectnb/siphot/howard/git/grating_synthesis'));

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

geometry = struct(  'thick_Si', thick_Si, ...
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

fprintf('Silicon Waveguide Thickness: %f\n', thick_Si);
fprintf('Silicon to Bottom Spacing: %f\n', space_Si_SiN);
fprintf('Bottom Nitride Thickness: %f\n', thick_bot_SiN);
fprintf('Bottom Nitride to Top Nitride Spacing: %f\n', space_SiN_SiN);
fprintf('Top Nitride Thickness: %f\n\n', thick_top_SiN);


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

% display Q for logging purposes
synth_obj

% design space fills for top bot version
fill_bots = 0.02:0.02:0.98;
fill_tops = 0.02:0.02:0.98;

% run design space generation, top bot version
tic;
synth_obj = synth_obj.generate_design_space(fill_bots, fill_tops);
toc;

fprintf('Design space generation sweep is done\n');
fprintf('Saving data...\n');

% make folder to save to
save_data_path = [  pwd filesep 'AIM'...
                    '__custom_thicknesses_' num2str(thick_Si) ...
                    '_' num2str(space_Si_SiN) ...
                    '_' num2str(thick_bot_SiN) ...
                    '_' num2str(space_SiN_SiN) ...
                    '_' num2str(thick_top_SiN) ...
                    '_nm__lam' num2str(lambda) ...
                    '__ang' num2str(optimal_angle) ...
                    '__dx_' num2str(synth_obj.discretization) ...
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

end


        
