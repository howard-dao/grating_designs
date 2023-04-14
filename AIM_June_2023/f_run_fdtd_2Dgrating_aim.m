function [ results ] = f_run_fdtd_2Dgrating_aim(synth_obj, geometry, ...
                                                single_or_dual, ...
                                                wg_interface_type, ...
                                                save_fname, OPTS)
% draw final design and simulate in fdtd
% for AIM stack

% inputs
%   synth_obj
%   geometry:
%       either 'a' or 'b'
%   single_or_dual
%       either 'single' or 'dual'
%   wg_interface_type
%       type: string
%       desc: waveguide interface type, currently supports 'body', 'poly', 
%       'body_and_poly', and 'nitride'
%   save_fname
%       type: string
%       desc: filepath and filename to save fsp file to
%   OPTS
%       type: struct
%       desc: optional settings
%               dxy - overrides discretization (units m)
%               AR_width - AR tooth width (m)
%               AR_offset - distance from left edge of AR tooth to grating
%       
% Outputs:
%   results
%       type: struct
%       desc: structure that contains results from the simulation + other
%               useful variables

% default si, custom sin, fully etched
[ box_thick, si_wg_thick, si_to_bot_sin_thick, bot_sin_thick, ...
        sin_to_sin_thick, top_sin_thick, ...
        n_SiO2, n_SiN, n_Si] = f_aim_layer_stack(   synth_obj.lambda * 1e-3, ...
                                                    geometry );

% convert units to m
box_thick           = 1e-6 * box_thick;
si_wg_thick         = 1e-6 * si_wg_thick;
si_to_bot_sin_thick = 1e-6 * si_to_bot_sin_thick;
bot_sin_thick       = 1e-6 * bot_sin_thick;
sin_to_sin_thick    = 1e-6 * sin_to_sin_thick;
top_sin_thick       = 1e-6 * top_sin_thick;

if strcmp( single_or_dual, 'dual' )
    synth_obj.synthesized_design.offset = 1e-9 * synth_obj.synthesized_design.offset;
end
synth_obj.synthesized_design.period = 1e-9 .* synth_obj.synthesized_design.period;

% mesh orders
mesh_ord_sio2       = 10;   % cladding
mesh_ord_si         = 5;    % silicon
mesh_ord_sin        = 8;    % nitride

% Start lumerical
% inputs
notes       = '';
lfs_name    = '2dfdtd.lsf';
[filepath,~,~] = fileparts(save_fname);

% make lumerical fdtd object
obj = c_lumericalFDTD(  'notes',            notes, ...
                        'filename',         lfs_name, ...
                        'file_directory',   filepath );

% open lumerical
obj = obj.appopen();

% grab wavelength
lambda0 = synth_obj.lambda * 1e-9;                         % units defined by synthobject

% grab synthesized design
synthesized_design = synth_obj.synthesized_design;

% get domain parameters (size, discretization, etc)
%     domain_size_y       = synthesized_design.y_coords(end) * synth_obj.units.scale;     % transverse size, in meters
domain_size_y       = 4e-6;
input_wg_length     = 6e-6;                                                         % length of input waveguide, m
output_gap_length   = 1e-6;                                                         % length of output space, m
domain_size_x       = input_wg_length + output_gap_length + synthesized_design.x_coords(end)*1e-9;
dxy                 = synth_obj.discretization * 1e-9;

if isfield(OPTS, 'dxy')
    dxy = OPTS.dxy;
end

mat_buff_size           = 5e-6; % area to extend outside of sim domain
monitor_displacement    = 2*dxy;

% add fdtd region
n_pml_layers    = 10;                                                          % number of pml discretizations
xmin            = 0;
xmax            = domain_size_x;
ymin            = 0;
ymax            = domain_size_y;
obj = obj.addfdtd(  'x min', xmin, ...
                    'x max', xmax, ...
                    'y min', ymin, ...
                    'y max', ymax, ...
                    'dimension', '2D', ...
                    'background index', synth_obj.background_index, ...
                    'mesh type', 'uniform', ...
                    'define x mesh by', 'maximum mesh step', ...
                    'define y mesh by', 'maximum mesh step', ...
                    'dx', dxy, ...
                    'dy', dxy, ...
                    'pml layers', n_pml_layers );
obj = obj.execute_commands();

% -------------------------
% Drawing structures

% draw layers

% draw input silicon layer
si_layer_miny = (ymin + ymax)/2 - si_wg_thick;
si_layer_maxy = si_layer_miny + si_wg_thick;
if strcmp( wg_interface_type, 'silicon')
    obj = obj.addrect(  'name',     'si layer', ...
                        'x min',    xmin - mat_buff_size, ...
                        'x max',    xmax + mat_buff_size, ...
                        'y min',    si_layer_miny, ...
                        'y max',    si_layer_maxy, ...
                        'z',        0, ...
                        'index',    n_Si, ...
                        'override mesh order from material database', true, ...
                        'mesh order', mesh_ord_si);
end

obj = obj.execute_commands();

% loop through synthesized design cells
cur_x           = xmin + input_wg_length;
bot_sin_miny    = si_layer_maxy + si_to_bot_sin_thick;
bot_sin_maxy    = bot_sin_miny + bot_sin_thick;
top_sin_miny    = bot_sin_maxy + sin_to_sin_thick;
top_sin_maxy    = top_sin_miny + top_sin_thick;

for ii = 1:length( synth.obj.synthesized_design.period )

    switch single_or_dual

        case 'single'
            % draw single level teeth

        case 'dual'
            % draw bi-level teeth

            % draw bottom tooth
            bot_tooth_len   = synth_obj.synthesized_design.period(ii) .* synth_obj.synthesized_design.bot_fill(ii);
            bot_tooth_minx  = cur_x;
            bot_tooth_maxx  = bot_tooth_minx + bot_tooth_len;
            obj = obj.addrect(  'name',     [ 'bottom_tooth_' num2str(ii) ], ...
                                'x min',    bot_tooth_minx, ...
                                'x max',    bot_tooth_maxx, ...
                                'y min',    bot_sin_miny, ...
                                'y max',    bot_sin_maxy, ...
                                'z',        0, ...
                                'index',    n_SiN, ...
                                'override mesh order from material database', true, ...
                                'mesh order', mesh_ord_sin );

            % draw top tooth
            top_tooth_len   = synth_obj.synthesized_design.period(ii) .* synth_obj.synthesized_design.top_fill(ii);
            top_tooth_minx  = bot_tooth_minx + synth_obj.synthesized_design.offset(ii);
            top_tooth_maxx  = top_tooth_minx + top_tooth_len;
            
            obj = obj.addrect(  'name',     [ 'top_tooth_' num2str(ii) ], ...
                                'x min',    top_tooth_minx, ...
                                'x max',    top_tooth_maxx, ...
                                'y min',    top_sin_miny, ...
                                'y max',    top_sin_maxy, ...
                                'z',        0, ...
                                'index',    n_SiN, ...
                                'override mesh order from material database', true, ...
                                'mesh order', mesh_ord_sin );
            % end case 'dual'

    end
    
    % update start x
    cur_x = cur_x + synth_obj.synthesized_design.period(ii);
end
obj = obj.execute_commands();

% -------------------------
% Adding sources and monitors

% add refractive index monitor
obj = obj.addindex( 'name', 'index_monitor', ...
                    'x min', 0, ...
                    'x max', domain_size_x, ...
                    'y min', 0, ...
                    'y max', domain_size_y );

% get index preview
[ obj, N ]      = obj.get_data_from_monitor( 'index_monitor', 'index preview' );
x               = N.index_monitor_index_preview.x;
y               = N.index_monitor_index_preview.y;
index_preview   = N.index_monitor_index_preview.index_x;                % isotropic (well except for pmls) so grabbing n_x
index_preview   = reshape( index_preview, length(x), length(y) ).';     % dimensions should be y, x

% add source
source_pos_x    = xmin + round(1e-6/dxy) * dxy;
obj             = obj.addmode( 'name', 'wg_mode_source', ...
                               'injection axis', 'x-axis', ...
                               'direction', 'Forward', ...
                               'mode selection', 'fundamental TM mode', ...
                               'x', source_pos_x, ...
                               'y min', 0, ...
                               'y max', domain_size_y, ...
                               'center wavelength', lambda0, ...
                               'wavelength span', 200e-9 );
% execute commands
obj = obj.execute_commands();

% add DFT field and power monitor across entire domain
n_freq_points = 1;
obj = obj.addpower( 'name', 'dft_full', ...
                    'simulation type', '2D Z-normal', ...
                    'monitor type', '2D Z-normal', ...
                    'x min', xmin, ...
                    'x max', xmax, ...
                    'y min', ymin, ...
                    'y max', ymax, ...
                    'override global monitor settings', true, ...
                    'use source limits', true, ...
                    'frequency points', n_freq_points );

% add DFT for transmission through upper region
n_freq_points = 40;
obj = obj.addpower( 'name', 'dft_up', ...
                    'simulation type', '2D Z-normal', ...
                    'monitor type', 'Linear X', ...
                    'x min', xmin, ...
                    'x max', xmax, ...
                    'y', ymax - monitor_displacement, ...
                    'override global monitor settings', true, ...
                    'use source limits', true, ...
                    'frequency points', n_freq_points ); 

% add DFT for transmission through bottom region
obj = obj.addpower( 'name', 'dft_down', ...
                    'simulation type', '2D Z-normal', ...
                    'monitor type', 'Linear X', ...
                    'x min', xmin, ...
                    'x max', xmax, ...
                    'y', ymin + monitor_displacement, ...
                    'override global monitor settings', true, ...
                    'use source limits', true, ...
                    'frequency points', n_freq_points ); 

% add DFT for transmission through end of domain
obj = obj.addpower( 'name', 'dft_throughput', ...
                    'simulation type', '2D Z-normal', ...
                    'monitor type', 'Linear Y', ...
                    'x', xmax - monitor_displacement, ...
                    'y min', ymin, ...
                    'y max', ymax, ...
                    'override global monitor settings', true, ...
                    'use source limits', true, ...
                    'frequency points', n_freq_points ); 

% add DFT for reflection
obj = obj.addpower( 'name', 'dft_reflection', ...
                    'simulation type', '2D Z-normal', ...
                    'monitor type', 'Linear Y', ...
                    'x', xmin + monitor_displacement, ...
                    'y min', ymin, ...
                    'y max', ymax, ...
                    'override global monitor settings', true, ...
                    'use source limits', true, ...
                    'frequency points', n_freq_points ); 
% execute commands
obj = obj.execute_commands();

% add mode overlap monitor for reflection
obj = obj.addmodeexpansion( ...
            'name', 'mode_expansion_monitor_r', ...
            'monitor type', '2D X-normal', ...
            'x', xmin + monitor_displacement, ...
            'y min', ymin, ...
            'y max', ymax, ...
            'mode selection', 'fundamental TM mode', ...
            'override global monitor settings', true, ...
            'use source limits', true, ...
            'frequency points', n_freq_points );
% add reflection monitor to the mode expansion monitor
obj = obj.setexpansion( 'output', 'dft_reflection' );
obj = obj.execute_commands();

% save the project
obj = obj.save( save_fname );
obj = obj.execute_commands();


% ------------------------------
% run simulation
obj = obj.run();
obj = obj.execute_commands();


% ------------------------------
% return results

% grab data from bottom monitor
[ obj, dft_down_T ] = obj.get_data_from_monitor( 'dft_down', 'T' );
lambda              = dft_down_T.dft_down_T.lambda;
frequencies         = dft_down_T.dft_down_T.f;
T_down              = abs( squeeze( dft_down_T.dft_down_T.T ) );
% Ez, Hx
[ obj, E_down ] = obj.get_data_from_monitor( 'dft_down', 'E' );
x_down          = E_down.dft_down_E.x;
Ez_down         = squeeze( E_down.dft_down_E.E(:,3,:) );    % dimensions field vs. x vs. freq
[ obj, H_down ] = obj.get_data_from_monitor( 'dft_down', 'H' );
Hx_down         = squeeze( H_down.dft_down_H.H(:,1,:) );    % dimensions field vs. x vs. freq

% grab data from upper monitor
[ obj, dft_up_T ]   = obj.get_data_from_monitor( 'dft_up', 'T' );
T_up                = abs( squeeze( dft_up_T.dft_up_T.T ) );
% Ez, Hx
[ obj, E_up ] = obj.get_data_from_monitor( 'dft_up', 'E' );
x_up          = E_up.dft_up_E.x;
Ez_up         = squeeze( E_up.dft_up_E.E(:,3,:) );    % dimensions field vs. x vs. freq
[ obj, H_up ] = obj.get_data_from_monitor( 'dft_up', 'H' );
Hx_up         = squeeze( H_up.dft_up_H.H(:,1,:) );    % dimensions field vs. x vs. freq

% grab the reflection transmission data
[ obj, dft_R ]  = obj.get_data_from_monitor( 'mode_expansion_monitor_r', 'expansion for output' );
R       = squeeze( dft_R.mode_expansion_monitor_r_expansion_for_output.T_total );
R_wg    = squeeze( dft_R.mode_expansion_monitor_r_expansion_for_output.T_backward );   % coupling into waveguide mode

% grab the throughport transmission data
[ obj, dft_thru ]   = obj.get_data_from_monitor( 'dft_throughput', 'T' );
T_thru              = abs( squeeze( dft_thru.dft_throughput_T.T ) );

% grab center wl (or closest to it)
[~, indx_center_wl]  = min( abs( lambda - lambda0 ) ); 
closest_lambda       = lambda( indx_center_wl );

% grab the correct data
if strcmp( synth_obj.coupling_direction, 'up' )
    Ez              = Ez_up;
    Hx              = Hx_up;
    Ez_center_freq  = Ez_up( :, indx_center_wl );
    Hx_center_freq  = Hx_up( :, indx_center_wl );
    T               = T_up;
    directivity     = 10*log10( T_up./T_down );                         % in dB
elseif strcmp( synth_obj.coupling_direction, 'down' )
    Ez              = Ez_down;
    Hx              = Hx_down;
    Ez_center_freq  = Ez_down( :, indx_center_wl );
    Hx_center_freq  = Hx_down( :, indx_center_wl );
    T               = T_down;
    directivity     = 10*log10( T_down./T_up );                         % in dB
end


% grab entire field data
[ obj, E_full_struct ]   = obj.get_data_from_monitor( 'dft_full', 'E' );
E_full      = E_full_struct.dft_full_E.E;
x_full      = E_full_struct.dft_full_E.x;
y_full      = E_full_struct.dft_full_E.y;
lambda_full = E_full_struct.dft_full_E.lambda;

% grab center wavelength (or closest to it)
[~, indx_center_wl_full]    = min( abs( lambda_full - lambda0 ) ); 
closest_lambda_full         = lambda_full( indx_center_wl_full );

Ez_all_center_freq  = E_full( :, 3, indx_center_wl_full );

% reshape the field
Ez_all_center_freq = squeeze(Ez_all_center_freq);
Ez_all_center_freq = reshape( Ez_all_center_freq, length(x_full), length(y_full) ).';


obj = obj.appclose();


% save results
results.lambda      = lambda;
results.frequencies = frequencies;

results.T_down  = T_down;
results.x_down  = x_down;
results.Ez_down = Ez_down;
results.Hx_down = Hx_down;

results.T_up    = T_up;
results.x_up    = x_up;
results.Ez_up   = Ez_up;
results.Hx_up   = Hx_up;

results.Ez_all_center_freq = Ez_all_center_freq;

results.R       = R;
results.R_wg    = R_wg;     % reflection into waveguide
results.dft_R   = dft_R;
results.T_thru  = T_thru;

results.directionality = directivity;

results.x = x;
results.y = y;
results.N = index_preview;

results.indx_center_wl = indx_center_wl;

results.single_or_dual      = single_or_dual;
%results.etch_depth          = etch_depth;
results.wg_interface_type   = wg_interface_type;
results.synth_obj           = synth_obj;

end