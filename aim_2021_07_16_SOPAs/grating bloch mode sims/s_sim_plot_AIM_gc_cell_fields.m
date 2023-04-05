% simulate aim cell with lower nitride gratings

clear; close all;

% dependencies
% laptop
addpath( genpath( [ 'C:\Users\howar\OneDrive\Desktop\Documents' ...
                    '\Boston University\Silicon Photonics\Grating Couplers' ...
                    '\grating_synthesis'] ) );
addpath( [  'C:\Users\howar\OneDrive\Desktop\Documents' ...
            '\Boston University\Silicon Photonics\Grating Couplers' ...
            '\grating_designs\aim_2021_07_16_SOPAs'] );

y_domain_size   = 4e3;
num_modes       = 3;
BC              = 0;    % 0 = PEC
pml_options     = [1, 100, 20, 2];
OPTS            = struct();

dxy     = 4;
lambda  = 1550;
k0      = 2*pi/lambda;

% neff_slab   = 2.842431701;
neff_slab   = 2.870152969854466;
n_clad      = 1.45;
theta       = 0; % deg
period      = 545;
period      = round(period/dxy)  * dxy;
guessk      = neff_slab*k0;

duty_cycles = [0, 0.275, 0.45, 0.816, 1];

figure;
for i_dc = 1:length(duty_cycles)
    fprintf('%d/%d Simulating for D_l=%g ...\n', i_dc, length(duty_cycles), duty_cycles(i_dc))

    % make grating cell then simulate
    GC = f_makeGratingCell_AIM( dxy, lambda, y_domain_size, period, 0, duty_cycles(i_dc), 0, OPTS );
    GC.numcells = 5;
    GC = GC.runSimulation( num_modes, BC, pml_options, k0, guessk);
    
    % plot index and field
    subplot(2,3,i_dc);
    x = 1e-3 * ( 0 : GC.dx : GC.domain_size(2)*GC.numcells - GC.dx );
    y = GC.y_coords*1e-3;
    N = repmat( GC.N, 1, GC.numcells );

    imagesc( x, y - y(end/2), real(GC.E_z) );
    title(['D_l=' num2str(duty_cycles(i_dc))]);
    colormap('redbluehilight');
    set( gca, 'ydir', 'normal' );
    set(gca, 'XTick', [], 'YTick', []);
    axis image;
    clim( 0.001 .* [-1,1] );
    ylim( [-1.5, 1.5] );
    % superimpose index contour
    hold on;
    contour( x, y - y(end/2), N, [2, 3.47], 'color', 'k', 'LineWidth', 1 );

end
%%

% include bilevel variant
GC = f_makeGratingCell_AIM( dxy, lambda, y_domain_size, period, 0.5, 0.74, 0.376147, OPTS );
GC.numcells = 5;
GC = GC.runSimulation( num_modes, BC, pml_options, k0, guessk);

subplot(2,3,6);
x = 1e-3 * ( 0 : GC.dx : GC.domain_size(2)*GC.numcells - GC.dx );
y = GC.y_coords*1e-3;
N = repmat( GC.N, 1, GC.numcells );

imagesc( x, y - y(end/2), real(GC.E_z) );
title(['D_l=' num2str(0.74) ', D_u=' num2str(0.5) ', O=' num2str(0.376)])
colormap('redbluehilight');
set( gca, 'ydir', 'normal' );
set(gca, 'XTick', [], 'YTick', []);
axis image;
clim( 0.001 .* [-1,1] );
ylim( [-1.5, 1.5] );
% superimpose index contour
hold on;
contour( x, y - y(end/2), N, [2, 3.47], 'color', 'k', 'LineWidth', 1 );



%%
decay_length = 1 / imag(GC.k) * 1e-6 / 2







    
    
    

