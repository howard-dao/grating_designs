function GC = f_makeGratingCell_AIM( dxy, lambda, y_domain_size, ...
                                     period, fill_top, fill_bot, offset_ratio, OPTS  )
% makes and returns a c_twoLevelGratingCell object
% for AIM process
% will eventually support customizable profile options
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
%       desc: wavelength to solve at, in units 'units'
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
%   OPTS
%       type: struct
%       desc: optional struct, going to hold parameters like "geometry"
%           fields: geometry - string
%
% outputs:
%   GC
%       type: c_twoLevelGratingCell object
%       desc: two level grating cell object

% if nargin < 7
%     OPTS = struct();
% end

% material indices
background_index    = 1.45;
n_Si                = 3.47;
n_SiN               = 2.0;

% geometry parameters (nm) (default parameters)
si_wg_thick         = 220;
sin_wg_thick        = 220;
si_to_sin_space     = 100;
sin_to_sin_space    = 100;

% set domain 
domain_size     = [ y_domain_size, period ];

% wrap offsets to range 0 to 1
offset_ratio = mod( offset_ratio, 1 );

% make 2 level grating cell
% GC = c_twoLevelGratingCell( 'discretization', dxy, ...
%                             'units', 'nm', ...
%                             'lambda', lambda, ...
%                             'domain_size', domain_size, ...
%                             'background_index', background_index, ...
%                             'numcells', 10 );
GC = c_twoLevelGratingCell( 'discretization', dxy, ...
                            'domain_size', domain_size, ...
                            'background_index', background_index, ...
                            'numcells', 10 );
                        
% halfway mark
domain_y_half   = round( (domain_size(1)/2) /GC.dy) * GC.dy;

if ~isfield( OPTS, 'geometry' )
    % default geometry

    % layer min heights
    si_miny     = domain_y_half - si_wg_thick; 
    si_maxy     = si_miny + si_wg_thick;
    botsin_miny = si_maxy + si_to_sin_space;
    botsin_maxy = botsin_miny + sin_wg_thick;
    topsin_miny = botsin_maxy + sin_to_sin_space;
    topsin_maxy = topsin_miny + sin_wg_thick;

    % add si layer
    GC = GC.addLayer(  si_miny, ...
                       si_wg_thick, ...
                       n_Si );   


    % draw two levels using two level builder function
    % the inputs are organized [ top level, bottom level ]
    wg_thick        = [ sin_wg_thick, sin_wg_thick ];
    wg_min_y        = [ topsin_miny, botsin_miny ];
    wgs_duty_cycles = [ fill_top, fill_bot ];
    wgs_offsets     = [ offset_ratio*period, 0 ];
    GC              = GC.twoLevelBuilder(   wg_min_y, ...
                                            wg_thick, ...
                                            [ n_SiN, n_SiN ], ...
                                            wgs_duty_cycles, ...
                                            wgs_offsets );

    % overriding min and max positions
    GC.wg_min_y = si_miny;
    GC.wg_max_y = topsin_maxy;
    
else
    
    switch( OPTS.geometry )
        
        case 'default si custom sin gratings full etch'
            % custom version that ovverides the sin thickness and spacings
            % but assumes fully etched sin layers (dual layers)
            % fields will be "si_to_bot_sin", "bot_sin_thick",
            % "sin_to_sin_thick", "top_sin_thick"
            
            
            % layer min heights
            si_miny     = domain_y_half - si_wg_thick; 
            si_maxy     = si_miny + si_wg_thick;
            botsin_miny = si_maxy + OPTS.si_to_bot_sin;
            botsin_maxy = botsin_miny + OPTS.bot_sin_thick;
            topsin_miny = botsin_maxy + OPTS.sin_to_sin_thick;
            topsin_maxy = topsin_miny + OPTS.top_sin_thick;

            % add si layer
            GC = GC.addLayer(  si_miny, ...
                               si_wg_thick, ...
                               n_Si );   


            % draw two levels using two level builder function
            % the inputs are organized [ top level, bottom level ]
            wg_thick        = [ OPTS.top_sin_thick, OPTS.bot_sin_thick ];
            wg_min_y        = [ topsin_miny, botsin_miny ];
            wgs_duty_cycles = [ fill_top, fill_bot ];
            wgs_offsets     = [ offset_ratio*period, 0 ];
            GC              = GC.twoLevelBuilder(   wg_min_y, ...
                                                    wg_thick, ...
                                                    [ n_SiN, n_SiN ], ...
                                                    wgs_duty_cycles, ...
                                                    wgs_offsets );

            % overriding min and max positions
            GC.wg_min_y = si_miny;
            GC.wg_max_y = topsin_maxy;
        
    end % end switch
    
end % end if/else
             
end




