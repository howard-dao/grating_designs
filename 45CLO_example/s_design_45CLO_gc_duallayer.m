% synthesize final designs for CLO

clear; close all;

% -----------------
% dependencies
addpath(genpath([   'C:\Users\howar\OneDrive\Desktop\Documents' ...
                    '\Boston University\Silicon Photonics\Grating Couplers\grating_synthesis']));
addpath(genpath(    'C:\Users\hdao\git\grating_synthesis'));
addpath(genpath(    '\projectnb\siphot\howard\git\grating_synthesis'));

addpath(genpath([   'C:\Users\howar\OneDrive\Desktop\Documents\' ...
                    'Boston University\Silicon Photonics\Grating Couplers\matlab_lumerical_api']))
addpath(genpath(    'C:\Users\hdao\git\matlab_lumerical_api'));
addpath(genpath(    '\projectnb\siphot\howard\git\matlab_lumerical_api'));

%addpath(genpath([   'C:\Users\howar\OneDrive\Desktop\Documents' ...
%                    '\Boston University\Silicon Photonics\Grating Couplers\my_matlab_utility_functions']));
%addpath(genpath('C:\Users\beezy\git\my_matlab_utility_functions'));

%addpath(genpath('C:\Program Files\Lumerical\v231\api\matlab'));

% -----------------

% synthesis object to load

% wl 1580, -15 deg, full etch, down
filepath = ['C:\Users\howar\OneDrive\Desktop\Documents' ...
            '\Boston University\Silicon Photonics\Grating Couplers\lam1580_ang-15_etch_full_dx_5_down'];
fill_top_override = 0.36;
fill_bot_override = 0.64;
fill_vs_top_bot_both = 'top'; % either be 'top', 'bottom', 'both'
etch_depth = 'full'; % 'full' or 'shallow'
MFD = 5e3; % mode field diameter, units are nm

% design options
apodized = true;

% fdtd options
run_fdtd = true;
fdtd_opts = struct('dxy', 10e-9);

f_design_45CLO_gc_duallayer(filepath, fill_top_override, fill_bot_override, ...
                            fill_vs_top_bot_both, etch_depth, MFD, apodized, ...
                            run_fdtd, fdtd_opts );




















