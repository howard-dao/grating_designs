% testing AIM custom grating cell

clear; close all;

% dependencies
addpath(genpath('C:\Users\howar\Siphot\Grating Couplers\grating_synthesis'));

dxy     = 10;
lambda  = 1550;
k0      = 2*pi/lambda;

neff_slab   = 2.842431701;
n_clad      = 1.45;
theta       = 15; % deg
period      = 2*pi/( k0 * ( neff_slab - n_clad*sin( theta*pi/180 ) ) );
period      = round(period/dxy)  * dxy;
y_domain_size = 4e3;

fill_top = 0.36;
fill_bot = 0.45;

offset_ratio = 0.4;

% custom layer thicknesses
thick_Si        = 220;
space_Si_SiN    = 100;
thick_bot_SiN   = 60;
space_SiN_SiN   = 0;
thick_top_SiN   = 160;

geometry = struct(  'thick_Si', thick_Si, ...
                    'space_Si_SiN', space_Si_SiN, ...
                    'thick_bot_SiN', thick_bot_SiN, ...
                    'space_SiN_SiN', space_SiN_SiN, ...
                    'thick_top_SiN', thick_top_SiN ...
                    );
           
GC = f_makeGratingCell_AIM_custom(  dxy, lambda, y_domain_size, period, ...
                                    fill_top, fill_bot, offset_ratio, geometry );

% plot index
GC.plotIndex();

% set grating solver settings
num_modes   = 15;
BC          = 0;    % 0 = PEC
pml_options = [1, 100, 20, 2]; 
OPTS        = struct();
guessk      = k0 * neff_slab;
% sim_opts    = struct('num_modes', num_modes, 'BC', BC, 'pml_options', pml_options);

% run simulation
GC = GC.runSimulation( num_modes, BC, pml_options, k0, guessk, OPTS );

% plot e field
GC.plot_E_field_gui();








