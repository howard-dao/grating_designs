function [n_SiO2, n_SiN, n_Si] = f_AIM_layer_indices(lambda_um)
% returns material index information for AIM process
%
% Inputs:
%   lambda
%       type: double, scalar
%       desc: wavelength [um]
%
% Outputs:
%   indices of refraction

% refractive indices
n_SiO2  = index_SiO2_fits(lambda_um);
n_SiN   = index_SiN(lambda_um);
n_Si    = index_Si_fits(lambda_um);

end

% -------------------------------------------------------------------------
% Auxiliary functions
% -------------------------------------------------------------------------

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

if nargin < 2
    temp = 298.15; % Default room temperature
end
lam = lam(:);
temp = temp(:);


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

if nargin < 2
    temp = 298.15; % Default room temperature
end
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
