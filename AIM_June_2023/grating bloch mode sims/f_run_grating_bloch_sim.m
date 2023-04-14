function [GC] = f_run_grating_bloch_sim(lambda, pitch, ...
                                   grating_type, duty_cycle, guessk, OPTS )
% runs the grating simulation
%
% inputs:
%   grating type
%       type: str
%       desc: either 'upper sin bars', 'lower sin bars', 'upper sin wg',
%       'lower sin wg', 'custom oxide thick sin bars'
%   OPTS
%       mode_to_overlap
%       custom_oxide_thickness
%
% for 'custom oxide thick sin bars', the oxide thickness between the
% silicon and the bottom-most sin layer is determined by
% custom_oxide_thickness (in nm). this only draws the bottom sin layer

if nargin < 6
    OPTS = struct();
end

% material indices
n_Si    = 3.48;
n_SiO2  = 1.45;
n_SiN   = 2.0;

% geometry parameters (nm)
si_wg_thick     = 220;
sin_wg_thick    = 220;
si_to_sin_space = 100;
BOX_thick       = 2000;

% grating cell settings
dxy                 = 10;   % nm
units               = 'nm';
domain_size         = [ 6000, pitch ];
background_index    = n_SiO2;
numcells            = 10;

% make single level grating cell
GC = c_gratingCell( 'discretization', dxy, ...
                    'units', units, ...
                    'lambda', lambda, ...
                    'domain_size', domain_size, ...
                    'background_index', background_index, ...
                    'numcells', numcells );
                
% draw the appropriate geometry
switch grating_type

    case 'upper sin bars'
        % si wg, upper sin bars
        
        % draw the silicon waveguide
        si_wg_bot = domain_size(1)/2 - si_wg_thick;
        GC = GC.addLayer(  si_wg_bot, ...
                           si_wg_thick, ...
                           n_Si );   
                       
        % draw upper nitride
        GC = GC.addRect( 0, ...
                         si_wg_bot + si_wg_thick + sin_wg_thick + 2*si_to_sin_space, ...
                         duty_cycle * pitch, ...
                         sin_wg_thick, ...
                         n_SiN );
                     
        % set wg boundaries
        GC.wg_min_y = si_wg_bot;
        GC.wg_max_y = si_wg_bot + si_wg_thick + si_to_sin_space + sin_wg_thick;
        
    case 'lower sin bars'
        
        % draw the silicon waveguide
        si_wg_bot = domain_size(1)/2 - si_wg_thick;
        GC = GC.addLayer(  si_wg_bot, ...
                           si_wg_thick, ...
                           n_Si );   
        
        GC = GC.addRect( 0, ...
                 si_wg_bot + si_wg_thick + si_to_sin_space, ...
                 duty_cycle * pitch, ...
                 sin_wg_thick, ...
                 n_SiN );
             
        % set wg boundaries
        GC.wg_min_y = si_wg_bot;
        GC.wg_max_y = si_wg_bot + si_wg_thick + si_to_sin_space + sin_wg_thick;
        
        
    case 'upper sin wg'
        
        % draw the sin waveguide
        si_wg_bot   = domain_size(1)/2 - si_wg_thick - sin_wg_thick - si_to_sin_space;
        sin_wg_bot  = si_wg_bot + si_wg_thick + sin_wg_thick + 2*si_to_sin_space;
        GC = GC.addLayer(  sin_wg_bot, ...
                           sin_wg_thick, ...
                           n_SiN );   
        
        % draw the si grating teeth
        GC = GC.addRect( 0, ...
                 si_wg_bot, ...
                 duty_cycle * pitch, ...
                 si_wg_thick, ...
                 n_Si );
             
        % set wg boundaries
        GC.wg_min_y = sin_wg_bot;
        GC.wg_max_y = sin_wg_bot + sin_wg_thick;
        
        
    case 'lower sin wg'
        
    case 'dual sin wg'
        
        % draw the sin waveguide
        si_wg_bot   = domain_size(1)/2 - si_wg_thick - sin_wg_thick - si_to_sin_space;
        sin_wg_bot  = si_wg_bot + si_wg_thick + si_to_sin_space;
        GC = GC.addLayer(  sin_wg_bot, ...
                           sin_wg_thick, ...
                           n_SiN );   
        GC = GC.addLayer(  sin_wg_bot + sin_wg_thick + si_to_sin_space, ...
                           sin_wg_thick, ...
                           n_SiN );  
        
        % draw the si grating teeth
        GC = GC.addRect( 0, ...
                 si_wg_bot, ...
                 duty_cycle * pitch, ...
                 si_wg_thick, ...
                 n_Si );
             
        % set wg boundaries
        GC.wg_min_y = sin_wg_bot;
        GC.wg_max_y = sin_wg_bot + sin_wg_thick*2 + si_to_sin_space;
    
    case 'custom oxide thick sin bars'
        
        % draw the silicon waveguide
        si_wg_bot = domain_size(1)/2 - si_wg_thick;
        GC = GC.addLayer(  si_wg_bot, ...
                           si_wg_thick, ...
                           n_Si );   
        
        % draw the sin grating tooth
        GC = GC.addRect( 0, ...
                 si_wg_bot + si_wg_thick + OPTS.custom_oxide_thickness, ...
                 duty_cycle * pitch, ...
                 sin_wg_thick, ...
                 n_SiN );
             
        % set wg boundaries
        GC.wg_min_y = si_wg_bot;
        GC.wg_max_y = si_wg_bot + si_wg_thick + OPTS.custom_oxide_thickness + sin_wg_thick;
        
    case 'unidir sin bars'
        
        % draw the silicon waveguide
        si_wg_bot = domain_size(1)/2 - si_wg_thick;
        GC = GC.addLayer(  si_wg_bot, ...
                           si_wg_thick, ...
                           n_Si );   
        
        % draw the sin grating tooth
        GC = GC.addRect( 0, ...
                 si_wg_bot + si_wg_thick + OPTS.custom_oxide_thickness, ...
                 duty_cycle * pitch, ...
                 sin_wg_thick, ...
                 n_SiN );
             
        % set wg boundaries
        GC.wg_min_y = si_wg_bot;
        GC.wg_max_y = si_wg_bot + si_wg_thick + OPTS.custom_oxide_thickness + sin_wg_thick;
        
end

% draw si substrate
si_sub_top      = si_wg_bot - BOX_thick;
GC              = GC.addLayer( 0, ...
                           si_sub_top, ...
                           n_Si );  

% debug, preview the index distribution
% GC.plotIndex();

% run simulation
num_modes   = 10;
BC          = 0;
pml_options = [1 200 20 2];
% guessk      = n_Si * (2 * pi / lambda);

% obj, num_modes, BC, pml_options, k0, guessk, OPTS );

if ~isfield(OPTS, 'mode_to_overlap')
    GC = GC.runSimulation( num_modes, BC, pml_options, 2*pi/lambda, guessk );
else
    GC = GC.runSimulation( num_modes, BC, pml_options, 2*pi/lambda, guessk, struct('mode_to_overlap',OPTS.mode_to_overlap) );
end

% debug, plot the field
% GC.plot_E_field_gui();

end



















