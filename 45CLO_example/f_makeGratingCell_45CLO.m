function GC = f_makeGratingCell_45CLO( dxy, lambda, y_domain_size, ...
                                         period, fill_top, fill_bot, offset_ratio, ...
                                         shallow_or_full )
% makes and returns a c_twoLevelGratingCell object
% with the 45CLO process parameters
%
% all units nm
%
% background index defaults to oxide
% 
% inputs:
%   dxy
%       type: double, scalar or 1x2 vector
%       desc: discretization along x and y
%             if scalar, then dx = dy = discretization
%             if vector, then discretization = [ dy dx ]
%   lambda
%       type: double, scalar
%       desc: wavelength
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
%   shallow_or_full
%       type: string
%       desc: 'shallow' for shallow etch, 'full' for full
%
% outputs:
%   GC
%       type: c_twoLevelGratingCell object
%       desc: two level grating cell object
%
% example:

% dependencies
% desktop
%addpath('C:\Users\bz\git\grating_synthesis\main');
%addpath('C:\Users\bz\git\grating_synthesis\chips\45CLO_2019_10\45CLO_functions');
% scc
%addpath('/projectnb/siphot/bz/git/grating_synthesis/main');
%addpath('/projectnb/siphot/bz/git/grating_synthesis/chips/45CLO_2019_10/45CLO_functions');

% Get process stack info, spatial units are in um
% switch units
%     case 'm'
%         lambda_um = lambda*1e6;
%     case 'mm'
%         lambda_um = lambda*1e3;
%     case 'um'
%         lambda_um = lambda;
%     case 'nm'
%         lambda_um = lambda*1e-3;
% end

% define layer thicknesses and refractive indices
[ t_box, t_bodysi, t_body_partial, t_polysi, t_nitride_liner, t_nitride_wg, ...
n_SiO2, n_SiN, n_cSi, n_pSi ] = f_clo_layer_stack( lambda*1e-3 );

% set domain 
domain_size = [ y_domain_size, period ];

% wrap offsets to range 0 to 1
offset_ratio = mod( offset_ratio, 1 );

% make 2 level grating cell
GC = c_twoLevelGratingCell( 'discretization', dxy, ...
                            'domain_size', domain_size, ...
                            'background_index', n_SiO2, ...
                            'numcells', 10 );

% convert units to nm
t_box           = 1e3 * t_box; 
t_bodysi        = 1e3 * t_bodysi;
t_body_partial  = 1e3 * t_body_partial;
t_polysi        = 1e3 * t_polysi;
t_nitride_liner = 1e3 * t_nitride_liner;

% define layer thicknesses
domain_y_half       = round( (domain_size(1)/2) /GC.dy) * GC.dy;
% t_air               = domain_y_half - t_box;

% layer min heights
body_miny = domain_y_half - t_bodysi;
poly_miny = domain_y_half;

% def. t_body_bottom as thickness of bottom half of body
% t_body_top as top half of body (thickness of bottom half of grating)
if strcmp( shallow_or_full, 'shallow' )
    t_body_bottom = t_body_partial;
elseif strcmp( shallow_or_full, 'full' )
    t_body_bottom = 0;
end
t_body_top = t_bodysi - t_body_bottom;

% draw layers
% GC = GC.addLayer( t_air, domain_size(1)-t_air, n_SiO2 );     % add in SiO2
% add body layer if shallow
if strcmp( shallow_or_full, 'shallow' )
    GC = GC.addLayer( body_miny, t_body_bottom, n_cSi );        % add in body layer, if shallow etch
end
GC = GC.addLayer( poly_miny, t_nitride_liner, n_SiN );          % add in SiN over body silicon
                        

% draw SiN rectangle, updated for offset on top tooth
top_layer_length = fill_top*period;
if top_layer_length > 0
                 
    % draw center tooth poly
    center_poly_tooth_minx  = offset_ratio*period;
    center_poly_tooth_maxx  = center_poly_tooth_minx + fill_top*period;
    nitride_minx            = center_poly_tooth_minx - t_nitride_liner;
    nitride_maxx            = center_poly_tooth_maxx + t_nitride_liner;
    nitride_miny            = poly_miny;
    nitride_maxy            = nitride_miny + t_polysi + t_nitride_liner;
    GC = GC.addRect( nitride_minx, ...
                     nitride_miny, ...
                     nitride_maxx - nitride_minx, ...
                     nitride_maxy - nitride_miny, ...
                     n_SiN );
                 
    % draw left tooth poly (assuming periodic)
    left_poly_tooth_minx  = center_poly_tooth_minx - period;
    left_poly_tooth_maxx  = left_poly_tooth_minx + fill_top*period;
    nitride_minx          = left_poly_tooth_minx - t_nitride_liner;
    nitride_maxx          = left_poly_tooth_maxx + t_nitride_liner;
    nitride_miny          = poly_miny;
    nitride_maxy          = nitride_miny + t_polysi + t_nitride_liner;
    GC = GC.addRect( nitride_minx, ...
                     nitride_miny, ...
                     nitride_maxx - nitride_minx, ...
                     nitride_maxy - nitride_miny, ...
                     n_SiN );
                 
    % draw right tooth poly
    right_poly_tooth_minx   = center_poly_tooth_minx + period;
    right_poly_tooth_maxx   = right_poly_tooth_minx + fill_top*period;
    nitride_minx            = right_poly_tooth_minx - t_nitride_liner;
    nitride_maxx            = right_poly_tooth_maxx + t_nitride_liner;
    nitride_miny            = poly_miny;
    nitride_maxy            = nitride_miny + t_polysi + t_nitride_liner;
    GC = GC.addRect( nitride_minx, ...
                     nitride_miny, ...
                     nitride_maxx - nitride_minx, ...
                     nitride_maxy - nitride_miny, ...
                     n_SiN );
                 
else
    % just a waveguide (body silicon)
    GC = GC.addLayer( poly_miny, t_nitride_liner, n_SiN );  % nitride layer   
end

% draw cell
% draw two levels using two level builder function
% the inputs are organized [ top level, bottom level ]
wg_thick        = [ t_polysi, t_body_top ];
wg_min_y        = [ poly_miny, body_miny + t_body_bottom ];
wgs_duty_cycles = [ fill_top, fill_bot ];
wgs_offsets     = [ offset_ratio*period, 0 ];
GC              = GC.twoLevelBuilder(   wg_min_y, wg_thick, [ n_pSi, n_cSi ], ...
                                        wgs_duty_cycles, wgs_offsets );
             
end


% -------------------------------------------------------------------------
% Auxiliary functions
% -------------------------------------------------------------------------



function n = index_IBM12SOI45nm_fits(Lum, material)
% Refractive index fits for IBM 12SOI 45nm process (used for EOS4 chip)
%
% author: unknown
%
% n = index_IBM12SOI45nm_fits(lam_um, material)
%
% Outputs:
% n        - [vector] Refractive index vector
%
% Inputs:
% lam_um   - [vector] wavelength in microns
% material - [string] 'polySi' - polysilicon gate
%                     'STI' - shallow trench isolation under waveguide
%                     'PSG' - oxide over the waveguide and Si3N4 liner

switch(material)
    case 'polySi'
        n = 4.337 - 1.247 * Lum + 0.6795 * Lum.^2 - 0.1305 * Lum.^3;
    case 'STI'
        n = 1.472 - 0.02556 * Lum + 0.005342 * Lum.^2;
    case 'PSG'
        n = 1.488 - 0.02998 * Lum + 0.00676 * Lum.^2;
    otherwise
        error([mfilename ': No valid material specified.']);
end

end         % end index_IBM12SOI45nm_fits()



function [nSi, nSimin, nSimax] = index_Si_fits(lam, temp)
% Silicon index vs. wavelength (1.5-1.6um) and temperature (85-920K)
% (simple fits of more complicated fits in source papers)
% author: M. Popovic, Apr 1, 2006
%
% Syntax:  [nSi nSimin nSimax] = index_Si_fits(lambda_um, temp_Kelvin)
%
% Input:    lambda_um      - [vector] wavelengths (accurate range 1.5-1.6um)
%           temp_Kelvin    - [vector] temperature (accurate range 85-920K)
%
% Output:   nSi            - nominal refractive index
%           nSimin, nSimax - 95% confidence interval of my fit of
%           literature fit curves.
%
% Sources of data:
% ----------------
% D.F. Edwards, "Silicon (Si)" in Palik, Handbook of Optical Consts, 1985.
%
% J.A. McCaulley et al., Phys. Rev. B 49, p7408 (Mar 1994): ‘Temperature
% dependence of the near-IR refractive index of Silicon, GaAs and InP’.

if(nargin < 2) temp = 298.15; end   % Default room temperature
lam = lam(:); temp = temp(:);

% Wavelength
A = -0.07805; dA = 0.0005;
B = 0.082; dB = 0.002;
x = lam - 1.55;
nn    = 3.476 + A*x + B*x.^2;               % Nominal n(lambda)
nnmin = 3.476 + (A-dA)*x + (B-dB)*x.^2;     % ..min and..
nnmax = 3.476 + (A+dA)*x + (B+dB)*x.^2;     % ..max of 95% conf interval

% Temperature
A = 0.01587; dA = 0.00001;
B = 0.002392; dB = 0.000005;
To = 298.15;  x = temp/To - 1;
nTovnTo = 1 + A*x + B*x.^2;
nTovnTomin = 1 + (A-dA)*x + (B-dB)*x.^2;
nTovnTomax = 1 + (A+dA)*x + (B+dB)*x.^2;

% Combine wavelength and temp, and set output
nSi = nn .* nTovnTo;  nSimin = nnmin .* nTovnTomin;  nSimax = nnmax .* nTovnTomax;

end         % end index_Si_fits()



function [nn, nnmin, nnmax] = index_SiO2_fits(lam, temp)
% Silica (SiO2) index vs. wavelength (x-xum) and temperature (x-xK)
% M. Popovic, Apr 28, 2006
%
% Syntax:  [n nmin nmax] = index_SiO2_fits(lambda_um, temp_Kelvin)
%
% Input:    lambda_um      - [vector] wavelengths (accurate range 1.5-1.6um)
%           temp_Kelvin    - [vector] temperature (accurate range 85-920K)
%
% Output:   n              - nominal refractive index
%           nmin, nmax     - 95% confidence interval of my fit of
%           literature fit curves.
%
% Sources of data:
% ----------------
% G. Ghosh et al., "Temperature-dependent Sellmeier coefficients and
% chromatic dispersions for some optical fiber glasses," J. Lightwave
% Technol. vol. 12, no. 8, Aug. 1994, p.1338.

if(nargin < 2) temp = 298.15; end   % Default room temperature
lam = lam(:); temp = temp(:);


TC = temp - 273.15;     % Celsius used in model
% Wavelength and temperature fit from above-cited paper
A = 6.90754e-6 * TC + 1.31552;
B = 2.35835e-5 * TC + 7.88404e-1;
C = 5.84758e-7 * TC + 1.10199e-2;
D = 5.48368e-7 * TC + 0.91316;
E = 100;
nn    = sqrt( A + B./(1 - C./lam.^2) + D./(1 - E./lam.^2) );  % Nominal n(lambda)
nnmin = NaN; %3.476 + (A-dA)*x + (B-dB)*x.^2;     % ..min and..
nnmax = NaN; %3.476 + (A+dA)*x + (B+dB)*x.^2;     % ..max of 95% conf interval

end         % end index_SiO2_fits()




function nn = index_SiN(lam)
% SiN refractive index fits to expt'al data
% http://www.filmetrics.com/refractive-index-database/Si3N4/Silicon-Nitride-SiN
%
% Syntax:  nn = index_SiN(lambda_um)
%
% Input:    lambda_um   [1-vector] wavelengths in microns; if blank, a plot is given.
%          
% Output:   nn          [1-vector] refractive index
%
% Cale Gentry, June 24, 2014

p = [-0.007287569871325   0.036401331486906  -0.079722168861525   2.052602158358897];

nn = p(1)*lam.^3+p(2)*lam.^2+p(3)*lam+p(4);

end





