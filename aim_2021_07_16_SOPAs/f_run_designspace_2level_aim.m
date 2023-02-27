function [] = f_run_designspace_2level_aim( lambda, optimal_angle, disc, ...
                                            coupling_direction, ...
                                            geometry, top_nitride_thick )
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
%   geometry
%       type: string
%       desc: layer stack settings: 'default' or 'custom'
%   top_nitride_thick
%       type: double, scalar
%       desc: thickness of top/second nitride layer [nm]

% dependencies
% desktop
addpath( genpath('C:\Users\hdao\git\grating_synthesis') );
% laptop
addpath( genpath(['C:\Users\howar\OneDrive\Desktop\Documents' ...
                    '\Boston University\Silicon Photonics' ...
                    '\Grating Couplers\grating_synthesis']))
% SCC
addpath( genpath('/projectnb/siphot/howard/git/grating_synthesis') );

% initial settings
units           = 'nm';
n_clad          = 1.45;
y_domain_size   = 4000;
data_notes      = [ 'lambda ', num2str(lambda), ...
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
fprintf('Designed for coupling into: %s\n', 'oxide index 1.45');
fprintf('Geometry: %s\n\n', geometry);

switch geometry
    case 'default'
        % default
        OPTS = struct();

        h_makeGratingCell = @(dxy, background_index, y_domain_size, ...
                             period, fill_top, fill_bot, offset_ratio) ...
                             f_makeGratingCell_AIM( dxy, ...
                                                    lambda, ...
                                                    y_domain_size, ...
                                                    period, ...
                                                    fill_top, ...
                                                    fill_bot, ...
                                                    offset_ratio , ...
                                                    OPTS );
    case 'custom'
        % default si, custom sin, fully etched
        si_to_bot_sin       = 100;
        bot_sin_thick       = 220;
        sin_to_sin_thick    = 50;
        top_sin_thick       = top_nitride_thick;

        OPTS = struct( 'geometry','default si custom sin gratings full etch', ...
                       'si_to_bot_sin', si_to_bot_sin, ...
                       'bot_sin_thick', bot_sin_thick, ...
                       'sin_to_sin_thick', sin_to_sin_thick, ...
                       'top_sin_thick', top_sin_thick );

        h_makeGratingCell = @(dxy, background_index, y_domain_size, ...
                             period, fill_top, fill_bot, offset_ratio) ...
                             f_makeGratingCell_AIM( dxy, ...
                                                    lambda, ...
                                                    y_domain_size, ...
                                                    period, ...
                                                    fill_top, ...
                                                    fill_bot, ...
                                                    offset_ratio, ...
                                                    OPTS  );
        fprintf('Silicon to Bottom Spacing: %f\n', si_to_bot_sin);
        fprintf('Bottom Nitride Thickness: %f\n', bot_sin_thick);
        fprintf('Bottom Nitride to Top Nitride Spacing: %f\n', sin_to_sin_thick);
        fprintf('Top Nitride Thickness: %f\n\n', top_sin_thick);
end

% make synthesis object
synth_obj = c_synthTwoLevelGrating(   'discretization',    disc, ...
                                      'units',             units,   ...
                                      'lambda',            lambda, ...
                                      'background_index',  n_clad,    ...
                                      'y_domain_size',     y_domain_size, ...
                                      'optimal_angle',     optimal_angle, ...
                                      'data_notes',        data_notes, ...
                                      'coupling_direction', coupling_direction, ...
                                      'h_makeGratingCell', h_makeGratingCell ...
                                      );

% display Q for logging purposes
synth_obj

% design space fills for top bot version
fill_bots = 0.06:0.02:0.96;
fill_tops = 0.06:0.02:0.96;

% run design space generation, top bot version
tic;
synth_obj = synth_obj.generate_design_space( fill_bots, fill_tops );
toc;

fprintf('Design space generation sweep is done\n');
fprintf('Saving data...\n');

% make folder to save to
switch geometry
    case 'default'
        save_data_path = [  pwd filesep 
                            'aimgc_' geometry ...
                            '_lambda' num2str(lambda) ...
                            '_optangle' num2str(optimal_angle) ...
                            '_dx_' num2str(synth_obj.discretization) ...
                            '_' coupling_direction ...
                            '_clad_' strrep(num2str(n_clad), '.', 'd') ];
    case 'custom'
        save_data_path = [  pwd filesep 
                            'aimgc_' geometry ...
                            '_top_nitride_thick' num2str(top_nitride_thick) ...
                            '_lambda' num2str(lambda) ...
                            '_optangle' num2str(optimal_angle) ...
                            '_dx_' num2str(synth_obj.discretization) ...
                            '_' coupling_direction ...
                            '_clad_' strrep(num2str(n_clad), '.', 'd') ];
mkdir( save_data_path );

% clear the GC from the data and save
save( [ save_data_path filesep 'synth_obj_aimgc' ], 'synth_obj' );

fprintf('Data saved\n');

% make folder to save plots to
save_plots_path = [ save_data_path filesep 'plots' ];
mkdir( save_plots_path );

fig_suffix = ''; %[ 'lambda' num2str(lambda) '_optangle' num2str(optimal_angle) '_etch_' shallow_or_full ];

f_save_design_space_plots( synth_obj, fig_suffix, save_plots_path );


end


        
