% simulates slab mode and picks wg thicknesses

clear; close all;

% dependencies
% laptop
addpath('C:\Users\howar\Siphot\grating_designs\slab_modesolver');

% desktop
addpath(genpath('C:\Users\hdao\Siphot\Grating Couplers\grating_designs\slab_modesolver'));

% wave parameters
lambda0     = 1550;              % um
k0          = 2*pi/lambda0;

% define material layers
n_clad      = 1.45;
n_core      = 3.477;
n           = [ n_clad, n_core, n_clad ];

% discretization
dx  = 10;

% define layer thicknesses
thick_clad      = 2e3;
thick_Si        = 220;
space_Si_SiN    = 100;
thick_bot_SiN   = 65;
space_SiN_SiN   = 100;

layer_t     = [ thick_clad, thick_Si, thick_clad ];

% choose # of modes to solve for
n_modes     = 5;

% pml options
% [ yes/no, length, strength, poly order ]
pml_opts = [ 0 0.2 10 2 ];

% guess k
guessk = 2*pi*n_core/lambda0;

% solve
tic;
[ field, k_fdfd, y, n_array ]  = f_slab_modesolver_bz( lambda0, n, layer_t, dx, guessk, n_modes, pml_opts );
toc;

% plot fundamental mode
figure;
plot( y, real(field(:,1) ) );
title('Fundamental Mode Field Slice');
xlabel('y Position [nm]'); ylabel('Field');
makeFigureNice();

% first figure out what the bottom overlap amount is
y = y - mean(y);

si_bot_sin_boty = thick_Si/2 + space_Si_SiN;
si_bot_sin_topy = si_bot_sin_boty + thick_bot_SiN;
total_bot_overlap = sum( field( y >= si_bot_sin_boty & y <= si_bot_sin_topy, 1 ) );

% plot the field within where the first nitride layer would be
figure;
plot( y( y >= si_bot_sin_boty & y <= si_bot_sin_topy ), field( y >= si_bot_sin_boty & y <= si_bot_sin_topy, 1 ) );
title('Field within Nitride Layer');
xlabel('y Position [nm]'); ylabel('Field');
makeFigureNice();

% debug
% total_bot_overlap_cumsum = 

% now get cumsum starting from where the top nitride should begin
si_top_sin_boty     = si_bot_sin_topy + space_SiN_SiN;
indx_topsin         = y >= si_top_sin_boty;
total_top_overlap_vs_topsin_thick = cumsum( field( indx_topsin, 1 ) );

figure;
plot( y( indx_topsin ) - min(y(indx_topsin)), total_top_overlap_vs_topsin_thick ); hold on;
plot( xlim, [ total_bot_overlap, total_bot_overlap ] );
xlabel('top SiN wg thick [nm]'); ylabel('Overlap');
legend('top wg overlap', 'target overlap');
title(['bot sin wg thick: ' num2str(thick_bot_SiN) ' nm, si to bot sin space: ' num2str(space_Si_SiN) ' nm']);
makeFigureNice();



























