function [ results ] = f_run_fdtd_2Dgrating_clo( synth_obj, etch_depth, ...
                                                single_or_dual, wg_interface_type, save_fname, OPTS )
% draw final design and simulate in fdtd
% for CLO stack
%
% inputs
%   synth_obj
%   etch_depth:
%       either 'full' or 'partial' or 'none'
%   single_or_dual
%       either 'single' or 'dual'
%   wg_interface_type
%       type: string
%       desc: waveguide interface type, currently supports 'body' and
%       'poly', 'body_and_poly', and 'nitride'
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
   

% % Choose whether single layer or dual layer grating
% single_or_dual = 'dual';  % either 'single' or 'dual'

if nargin < 6
    OPTS = struct();
end

% get process info
[ t_box, t_bodysi, t_body_partial, t_polysi, t_nitride_liner, t_nitride_wg, ...
           n_SiO2, n_SiN, n_cSi, n_pSi ] = f_clo_layer_stack( synth_obj.lambda * 1e-3 );

% convert units to m
t_box           = 1e-6 * t_box; 
t_bodysi        = 1e-6 * t_bodysi;
t_body_partial  = 1e-6 * t_body_partial;
t_polysi        = 1e-6 * t_polysi;
t_nitride_liner = 1e-6 * t_nitride_liner;
t_nitride_wg    = 1e-6 * t_nitride_wg;

if strcmp( single_or_dual, 'dual' )
    synth_obj.synthesized_design.offset = 1e-9 * synth_obj.synthesized_design.offset;
end
synth_obj.synthesized_design.period = 1e-9 .* synth_obj.synthesized_design.period;
       
% mesh orders
mesh_ord_sio2       = 10;   % cladding
mesh_ord_si_body    = 5;    % silicon body layer
mesh_ord_sin        = 8;    % nitride
mesh_ord_si_poly    = 7;    % poly
mesh_ord_sietch     = 4;    % si etch
mesh_ord_sitooth    = 3;    % silicon tooth
mesh_ord_nitride_wg = 2;    % nitride waveguide

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

% draw silicon body layer
si_layer_miny = (ymin + ymax)/2 - t_bodysi;
si_layer_maxy = si_layer_miny + t_bodysi;
if strcmp( wg_interface_type, 'body' ) || strcmp( wg_interface_type, 'body_and_poly' )
    % Draw body layer
    obj = obj.addrect( 'name', 'si body layer', ...
                       'x min', xmin - mat_buff_size, ...
                       'x max', xmax + mat_buff_size, ...
                       'y min', si_layer_miny, ...
                       'y max', si_layer_maxy, ...
                       'z',     0, ...
                       'index', n_cSi, ...
                       'override mesh order from material database', true, ...
                       'mesh order', mesh_ord_si_body );
end

% draw nitride liner layer
sin_layer_miny = si_layer_maxy;
sin_layer_maxy = sin_layer_miny + t_nitride_liner;
obj = obj.addrect( 'name', 'sin layer', ...
                   'x min', xmin - mat_buff_size, ...
                   'x max', xmax + mat_buff_size, ...
                   'y min', sin_layer_miny, ...
                   'y max', sin_layer_maxy, ...
                   'z',     0, ...
                   'index', n_SiN, ...
                   'override mesh order from material database', true, ...
                   'mesh order', mesh_ord_sin );
obj = obj.execute_commands();

% draw poly layer, if waveguide interface type is poly
if strcmp( wg_interface_type, 'poly' )
    poly_layer_miny = si_layer_maxy;
    poly_layer_maxy = poly_layer_miny + t_polysi;
    obj = obj.addrect( 'name', 'si poly layer', ...
                       'x min', xmin, ...
                       'x max', xmax, ...
                       'y min', poly_layer_miny, ...
                       'y max', poly_layer_maxy, ...
                       'z',     0, ...
                       'index', n_pSi, ...
                       'override mesh order from material database', true, ...
                       'mesh order', mesh_ord_si_poly );
    % draw nitride liner layer
    sinp_layer_miny = poly_layer_miny;
    sinp_layer_maxy = poly_layer_maxy + t_nitride_liner;
    obj = obj.addrect( 'name', 'sin layer on poly', ...
                       'x min', xmin, ...
                       'x max', xmax, ...
                       'y min', sinp_layer_miny, ...
                       'y max', sinp_layer_maxy, ...
                       'z',     0, ...
                       'index', n_SiN, ...
                       'override mesh order from material database', true, ...
                       'mesh order', mesh_ord_sin );
    obj = obj.execute_commands();
end

% draw nitride waveguide layer, if waveguide interface is nitride
nitride_wg_miny = sin_layer_maxy;
nitride_wg_maxy = nitride_wg_miny + t_nitride_wg;
if strcmp( wg_interface_type, 'nitride' )
    % Draw nitride layer
    obj = obj.addrect( 'name', 'nitride layer', ...
                       'x min', xmin, ...
                       'x max', xmax, ...
                       'y min', nitride_wg_miny, ...
                       'y max', nitride_wg_maxy, ...
                       'z',     0, ...
                       'index', n_SiN, ...
                       'override mesh order from material database', true, ...
                       'mesh order', mesh_ord_nitride_wg );
    obj = obj.execute_commands();
end
    
% etch waveguide
etch_minx  = xmin + input_wg_length;
etch_maxx  = xmax + mat_buff_size;
if ( strcmp( etch_depth, 'partial' ) || strcmp( etch_depth, 'shallow' ) )
    % draw partial etch
    petch_miny  = si_layer_miny + t_body_partial;
    petch_maxy  = si_layer_maxy;
    
    obj = obj.addrect( 'name', 'partial etch layer', ...
                       'x min', etch_minx, ...
                       'x max', etch_maxx, ...
                       'y min', petch_miny, ...
                       'y max', petch_maxy, ...
                       'z',     0, ...
                       'index', n_SiO2, ...
                       'override mesh order from material database', true, ...
                       'mesh order', mesh_ord_sietch );
    obj = obj.execute_commands();
elseif strcmp( etch_depth, 'full' )
    % draw full etch
    etch_miny  = si_layer_miny;
    etch_maxy  = si_layer_maxy;
    obj = obj.addrect( 'name', 'full etch layer', ...
                       'x min', etch_minx, ...
                       'x max', etch_maxx, ...
                       'y min', etch_miny, ...
                       'y max', etch_maxy, ...
                       'z',     0, ...
                       'index', n_SiO2, ...
                       'override mesh order from material database', true, ...
                       'mesh order', mesh_ord_sietch );
    obj = obj.execute_commands();
end
    
% add AR tooth
if isfield(OPTS, 'AR_width')
   
    % currently only supporting RX full etch
    AR_minx  = input_wg_length - OPTS.AR_offset;
    AR_maxx  = AR_minx + OPTS.AR_width;
    obj = obj.addrect( 'name', 'AR tooth', ...
                       'x min', AR_minx, ...
                       'x max', AR_maxx, ...
                       'y min', etch_miny, ...
                       'y max', etch_maxy, ...
                       'z',     0, ...
                       'index', n_SiO2, ...
                       'override mesh order from material database', true, ...
                       'mesh order', mesh_ord_sietch );
    obj = obj.execute_commands();
    
end
    
% % add thick input waveguide, if that was the design
% if strcmp( synthesized_design.input_wg_type, 'full' ) && ~strcmp( wg_interface_type, 'nitride' )
% 
%     % add poly si layer
%     polysi_miny = si_layer_maxy;
%     polysi_maxy = polysi_miny + t_polysi;
%     obj = obj.addrect( 'name', 'polysi layer', ...
%                    'x min', xmin - mat_buff_size, ...
%                    'x max', xmin + input_wg_length, ...
%                    'y min', polysi_miny, ...
%                    'y max', polysi_maxy, ...
%                    'z',     0, ...
%                    'index', n_pSi, ...
%                    'override mesh order from material database', true, ...
%                    'mesh order', mesh_ord_si_poly );
% 
%     % add poly nitride liner
% %         polysin_maxx = input_wg_length + t_nitride_liner;
%     polysin_miny = si_layer_maxy;
%     polysin_maxy = polysi_maxy + t_nitride_liner;
%     obj = obj.addrect( 'name', 'polysi nitride layer', ...
%                    'x min', xmin - mat_buff_size, ...
%                    'x max', xmin + input_wg_length + t_nitride_liner, ...
%                    'y min', polysin_miny, ...
%                    'y max', polysin_maxy, ...
%                    'z',     0, ...
%                    'index', n_SiN, ...
%                    'override mesh order from material database', true, ...
%                    'mesh order', mesh_ord_sin );
% 
% end

% set bottom tooth thickness
if ( strcmp( etch_depth, 'partial' ) || strcmp( etch_depth, 'shallow' ) )
    % partial etch
    body_tooth_miny = petch_miny;
    body_tooth_maxy = petch_maxy;
elseif strcmp( etch_depth, 'full' )
    % full etch
    body_tooth_miny = si_layer_miny;
    body_tooth_maxy = si_layer_maxy;
end

% set top tooth thickness (poly)
si_poly_miny    = si_layer_maxy;
si_poly_maxy    = si_poly_miny + t_polysi;

% set nitride thickness
polysin_miny = si_poly_miny;
polysin_maxy = si_poly_maxy + t_nitride_liner;


% loop through synthesized design cells
cur_x   = xmin + input_wg_length;
for ii = 1:length( synth_obj.synthesized_design.period )

    switch single_or_dual

        case 'single'
            % draw single level teeth

            if strcmp( wg_interface_type, 'body' ) && strcmp( etch_depth, 'none' )
                % poly teeth, input waveguide was body and etch depth
                % was none
                si_poly_minx    = cur_x;
                si_poly_maxx    = si_poly_minx + synth_obj.synthesized_design.period(ii) .* synth_obj.synthesized_design.fill(ii);
                si_poly_miny    = si_layer_maxy;
                si_poly_maxy    = si_poly_miny + t_polysi;
                obj = obj.addrect( 'name', [ 'si_poly_' num2str(ii)], ...
                   'x min', si_poly_minx, ...
                   'x max', si_poly_maxx, ...
                   'y min', si_poly_miny, ...
                   'y max', si_poly_maxy, ...
                   'z',     0, ...
                   'index', n_pSi, ...
                   'override mesh order from material database', true, ...
                   'mesh order', mesh_ord_si_poly );
                % nitride on top of poly
                polysin_minx = si_poly_minx - t_nitride_liner;
                polysin_maxx = si_poly_maxx + t_nitride_liner;
                polysin_miny = si_poly_miny;
                polysin_maxy = si_poly_maxy + t_nitride_liner;
                obj = obj.addrect( 'name', [ 'polysin_' num2str(ii)], ...
                   'x min', polysin_minx, ...
                   'x max', polysin_maxx, ...
                   'y min', polysin_miny, ...
                   'y max', polysin_maxy, ...
                   'z',     0, ...
                   'index', n_SiN, ...
                   'override mesh order from material database', true, ...
                   'mesh order', mesh_ord_sin );

            elseif ( strcmp( wg_interface_type, 'poly' ) && strcmp( etch_depth, 'full' ) ) || strcmp( wg_interface_type, 'nitride' )
                % body teeth, input waveguide was poly
                % was full
                % OR nitride interfacing waveguide with body teeth
                si_body_minx    = cur_x;
                si_body_maxx    = si_body_minx + synth_obj.synthesized_design.period(ii) .* synth_obj.synthesized_design.fill(ii);
                si_body_miny    = si_layer_miny;
                si_body_maxy    = si_body_miny + t_bodysi;
                obj = obj.addrect( 'name', [ 'si_body_' num2str(ii)], ...
                   'x min', si_body_minx, ...
                   'x max', si_body_maxx, ...
                   'y min', si_body_miny, ...
                   'y max', si_body_maxy, ...
                   'z',     0, ...
                   'index', n_cSi, ...
                   'override mesh order from material database', true, ...
                   'mesh order', mesh_ord_sitooth );
               
            elseif ( strcmp( wg_interface_type, 'body' ) && ( strcmp( etch_depth, 'partial' ) || strcmp( etch_depth, 'full' ) ) )
                % single layer shallow etch
                % actually this would work for full etch too
                obj = obj.addrect( 'name', [ 'si_tooth_' num2str(ii)], ...
                   'x min', cur_x, ...
                   'x max', cur_x + synth_obj.synthesized_design.period(ii) .* synth_obj.synthesized_design.fill(ii), ...
                   'y min', si_layer_miny, ...
                   'y max', si_layer_miny + t_bodysi, ...
                   'z',     0, ...
                   'index', n_cSi, ...
                   'override mesh order from material database', true, ...
                   'mesh order', mesh_ord_sitooth );
               
            end     % end if strcmp( wg_interface_type, 'body' ) && strcmp( etch_depth, 'none' )

        case 'dual'
            % draw bi level teeth
            
            % position
            body_tooth_minx = cur_x;
            body_tooth_maxx = body_tooth_minx + synth_obj.synthesized_design.period(ii) .* synth_obj.synthesized_design.bot_fill(ii);

            % draw
            obj = obj.addrect( 'name', [ 'bodytooth_' num2str(ii)], ...
               'x min', body_tooth_minx, ...
               'x max', body_tooth_maxx, ...
               'y min', body_tooth_miny, ...
               'y max', body_tooth_maxy, ...
               'z',     0, ...
               'index', n_cSi, ...
               'override mesh order from material database', true, ...
               'mesh order', mesh_ord_sitooth );

            % draw the top tooth
            % poly silicon layer
            si_poly_len     = synth_obj.synthesized_design.period(ii) .* synth_obj.synthesized_design.top_fill(ii);
            si_poly_minx    = body_tooth_minx + synth_obj.synthesized_design.offset(ii);
            si_poly_maxx    = si_poly_minx + si_poly_len;
            

             % nitride on top of poly
            polysin_minx = si_poly_minx - t_nitride_liner;
            polysin_maxx = si_poly_maxx + t_nitride_liner;

            obj = obj.addrect( 'name', [ 'si_poly_' num2str(ii)], ...
               'x min', si_poly_minx, ...
               'x max', si_poly_maxx, ...
               'y min', si_poly_miny, ...
               'y max', si_poly_maxy, ...
               'z',     0, ...
               'index', n_pSi, ...
               'override mesh order from material database', true, ...
               'mesh order', mesh_ord_si_poly );

            obj = obj.addrect( 'name', [ 'polysin_' num2str(ii)], ...
               'x min', polysin_minx, ...
               'x max', polysin_maxx, ...
               'y min', polysin_miny, ...
               'y max', polysin_maxy, ...
               'z',     0, ...
               'index', n_SiN, ...
               'override mesh order from material database', true, ...
               'mesh order', mesh_ord_sin );

       % end case 'dual'

    end     % end switch

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
results.etch_depth          = etch_depth;
results.wg_interface_type   = wg_interface_type;
results.synth_obj           = synth_obj;

end

% -------------------------------------------------------------------------
% Aux functions
% -------------------------------------------------------------------------

















