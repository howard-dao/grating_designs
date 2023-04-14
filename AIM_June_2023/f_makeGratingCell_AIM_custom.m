function GC = f_makeGratingCell_AIM_custom( dxy, lambda, y_domain_size, period, ...
                                            fill_top, fill_bot, offset_ratio, geometry )
% makes and returns a c_twoLevelGratingCell object
% for AIM process with custom layer stack
% 
% authors: howard
% 
% units are in nm, specifically for this method
% 
% inputs:
%   dxy
%       type: double, scalar or 1x2 vector
%       desc: discretization along x and y, in units of 'units'
%             if scalar, then dx = dy = discretization
%             if vector, then discretization = [ dy dx ]
%   lambda
%       type: double, scalar
%       desc: wavelength to solve at, in nm
%   y_domain_size
%       type: double, scalar
%       desc: vertical/transverse domain size
%   period
%       type: double, scalar
%       desc: period of the grating cell, in units defined by synth_obj.units
%   fill_top
%       type: double, scalar
%       desc: ratio of top layer to period
%   fill_bot
%       type: double, scalar
%       desc: ratio of bottom layer to bottom layer
%   offset_ratio
%       type: double, scalar
%       desc: ratio of TOP layer offset to period
%   geometry
%       type: struct
%       desc: parameters for layer thicknesses
% 
% outputs:
%   GC
%       type: c_twoLevelGratingCell object
%       desc: two level grating cell object

% material indices
% n_SiO2  = 1.445;
% n_Si    = 3.476;
% n_SiN   = 2.0;
[n_SiO2, n_SiN, n_Si] = f_AIM_layer_indices(1e-3 * lambda);

% geometry parameters [nm]
thick_Si        = geometry.thick_Si;
space_Si_SiN    = geometry.space_Si_SiN;
thick_bot_SiN   = geometry.thick_bot_SiN;
space_SiN_SiN   = geometry.space_SiN_SiN;
thick_top_SiN   = geometry.thick_top_SiN;

% set domain
domain_size     = [y_domain_size, period];

% wrap offsets to range 0 to 1
offset_ratio = mod( offset_ratio, 1 );

% make 2 level grating cell
GC = c_twoLevelGratingCell( 'discretization', dxy, ...
                            'domain_size', domain_size, ...
                            'background_index', n_SiO2, ...
                            'numcells', 10 ...
                            );

% halfway mark
domain_y_half = round( (domain_size(1)/2) / GC.dy) * GC.dy;

% layer min heights
si_miny     = domain_y_half - thick_Si;
si_maxy     = si_miny + thick_Si;
botsin_miny = si_maxy + space_Si_SiN;
botsin_maxy = botsin_miny + thick_bot_SiN;
topsin_miny = botsin_maxy + space_SiN_SiN;
topsin_maxy = topsin_miny + thick_top_SiN;

% add si layer
GC = GC.addLayer(  si_miny, ...
                   thick_Si, ...
                   n_Si );   

% draw two levels using two level builder function
% the inputs are organized [ top level, bottom level ]
wg_thick        = [thick_top_SiN, thick_bot_SiN];
wg_min_y        = [topsin_miny, botsin_miny];
wgs_duty_cycles = [fill_top, fill_bot];
wgs_offsets     = [offset_ratio * period, 0];
GC              = GC.twoLevelBuilder(   wg_min_y, ...
                                        wg_thick, ...
                                        [ n_SiN, n_SiN ], ...
                                        wgs_duty_cycles, ...
                                        wgs_offsets ...
                                        );

% overriding min and max positions
GC.wg_min_y = si_miny;
GC.wg_max_y = topsin_maxy;




