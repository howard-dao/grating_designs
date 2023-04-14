% synthesize final designs for AIM

clear; close all;

% -----------------
% dependencies
% laptop
addpath(genpath('C:\Users\howar\Siphot\Grating Couplers\grating_synthesis'));
addpath(genpath('C:\Users\howar\Siphot\Grating Couplers\matlab_lumerical_api'));

% desktop
addpath(genpath('C:\Users\hdao\git\grating_synthesis'));
addpath(genpath('C:\Users\hdao\git\matlab_lumerical_api'));

% SCC
addpath(genpath('\projectnb\siphot\howard\git\grating_synthesis'));
addpath(genpath('\projectnb\siphot\howard\git\matlab_lumerical_api'));

% lumerical
addpath(genpath('C:\Program Files\Lumerical\v231\api\matlab'));

% -----------------

% synthesis object to load

% wl 1550, +15 deg, up, custom
filepath = ['C:\Users\howar\OneDrive\Desktop\Documents' ...
            '\Boston University\Silicon Photonics\Grating Couplers' ...
            '\grating_designs\aim_2021_07_16_SOPAs\bash_scripts\' ...
            '2023_02_17_17_04_37_aimgc_geometry_b_lambda1550_optangle15_dx_5_up_clad_1d45'];

fill_top_override = 0.35;
fill_bot_override = 0.50;
fill_vs_top_bot_both = 'bottom';    % either be 'top', 'bottom', 'both'
geometry = 'b';                     % either 'a' or 'b' for now
MFD = 5e3;                          % mode field diameter, units are nm

% design options
apodized = true;

% fdtd options
run_fdtd = false;
fdtd_opts = struct('dxy', 10e-9);

f_design_AIM_gc_duallayer(  filepath, fill_top_override, fill_bot_override, ...
                            fill_vs_top_bot_both, geometry, MFD, apodized, ...
                            run_fdtd, fdtd_opts );

