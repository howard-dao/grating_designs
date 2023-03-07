% simulates slab mode and picks wg thicknesses

clear; close all;

addpath([   'C:\Users\howar\OneDrive\Desktop\Documents' ...
            '\Boston University\Silicon Photonics\Grating Couplers' ...
            '\grating_designs\aim_2021_07_16_SOPAs\slab_modesolver']);
addpath( (genpath( [ 'C:\Users\howar\OneDrive\Desktop\Documents' ...
                     '\Boston University\Silicon Photonics\Grating Couplers' ...
                     '\grating_designs\aim_2021_07_16_SOPAs'])));

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
clad_d      = 2e3;
si_d        = 220;
layer_t     = [ clad_d, si_d, clad_d ];

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
xlabel('y Position (\mum)'); ylabel('Field');
makeFigureNice();

% now lets pick where the nitride should be
% fields will be "si_to_bot_sin", "bot_sin_thick",
        % "sin_to_sin_thick", "top_sin_thick"
sin_index           = 2.0;
si_to_bot_sin       = 100;
bot_sin_thick       = 220;
sin_to_sin_space    = 100;

% first figure out what the bottom overlap amount is
y = y - mean(y);

si_bot_sin_boty = si_d/2 + si_to_bot_sin;
si_bot_sin_topy = si_bot_sin_boty + bot_sin_thick;
total_bot_overlap = sum( field( y >= si_bot_sin_boty & y <= si_bot_sin_topy, 1 ) );

% plot the field within where the first nitride layer would be
figure;
plot( y( y >= si_bot_sin_boty & y <= si_bot_sin_topy ), field( y >= si_bot_sin_boty & y <= si_bot_sin_topy, 1 ) );
title('Field within Nitride Layer');
xlabel('y Position (\mum)'); ylabel('Field');
makeFigureNice();

% debug
% total_bot_overlap_cumsum = 

% now get cumsum starting from where the top nitride should begin
si_top_sin_boty     = si_bot_sin_topy + sin_to_sin_space;
indx_topsin         = y >= si_top_sin_boty;
total_top_overlap_vs_topsin_thick = cumsum( field( indx_topsin, 1 ) );

figure;
plot( y( indx_topsin ) - min(y(indx_topsin)), total_top_overlap_vs_topsin_thick ); hold on;
plot( xlim, [ total_bot_overlap, total_bot_overlap ] );
xlabel('top SiN wg thick (\mum)'); ylabel('Overlap');
legend('top wg overlap', 'target overlap');
title(['bot sin wg thick: ' num2str(bot_sin_thick) ' nm, si to bot sin space: ' num2str(si_to_bot_sin) ' nm']);
makeFigureNice();



























