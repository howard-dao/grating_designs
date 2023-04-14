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
dx  = 5;

% define layer thicknesses
thick_clad      = 2e3;
thick_Si        = 220;
space_Si_SiN    = 100;
thick_bot_SiN   = 0 : dx : 300;
% space_SiN_SiN   = 100;

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


% first figure out what the bottom overlap amount is
y = y - mean(y);

si_bot_sin_boty = thick_Si/2 + space_Si_SiN;

si_bot_sin_topy = si_bot_sin_boty + thick_bot_SiN;

total_bot_overlap = zeros(size(si_bot_sin_topy));
for i_bot_thick = 1:length(thick_bot_SiN)
    total_bot_overlap(i_bot_thick) = sum(field( y >= si_bot_sin_boty & ...
                                                y <= si_bot_sin_topy(i_bot_thick), ...
                                                1 ) );
end

% plot bottom nitride overlap with unperturbed field as a function of
% bottom nitride thickness
figure;
plot(thick_bot_SiN, total_bot_overlap, 'Color', 'b', 'LineWidth', 1.5); hold on;
plot([65, 65], [0, 1], 'Color', 'k', 'LineWidth', 2); hold on;
plot([220, 220], [0, 1], 'Color', 'k', 'LineWidth', 2);
title('Bottom SiN Field Overlap with Unperturbed Field');
xlabel('Bottom SiN Thickness [nm]');
ylabel('Total Bottom Overlap');
grid on;






















