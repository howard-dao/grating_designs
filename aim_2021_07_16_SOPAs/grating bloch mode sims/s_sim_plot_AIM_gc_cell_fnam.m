% simulate aim cell for DII

clear; close all;

% dependencies
% laptop
addpath( genpath( [ 'C:\Users\howar\OneDrive\Desktop\Documents' ...
                    '\Boston University\Silicon Photonics\Grating Couplers' ...
                    '\grating_synthesis'] ) );

y_domain_size       = 4e3;
num_modes   = 5;
BC          = 0;    % 0 = PEC
pml_options = [1, 100, 20, 2];
OPTS        = struct();

% simulate and plot for a bottom nitride grating of similar duty cycle
duty_cycle = 0.78;

dxy     = 10;
lambda  = 1550;
k0      = 2*pi/lambda;

neff_slab   = 2.842431701;
n_clad      = 1.45;
theta       = 15; % deg
period      = 2*pi/( k0 * ( neff_slab - n_clad*sin( theta*pi/180 ) ) );
period      = round(period/dxy)  * dxy;
guessk      = neff_slab*k0;

GC_botsin = f_makeGratingCell_AIM(  dxy, lambda, y_domain_size, ...
                                    period, 0, duty_cycle, 0, OPTS );
GC_botsin.numcells = 5;
                           
% plot index
GC_botsin.plotIndex();

% run simulation
tic;
GC_botsin = GC_botsin.runSimulation( num_modes, BC, pml_options, k0, guessk, OPTS );
toc;

% plot e field
% GC.plot_E_field_gui();

% plot field with waveguide overlay
x = 1e-3 * ( 0 : GC_botsin.dx : GC_botsin.domain_size(2)*GC_botsin.numcells - GC_botsin.dx );
y = GC_botsin.y_coords*1e-3;
N = repmat( GC_botsin.N, 1, GC_botsin.numcells );

figure;
imagesc( x, y - y(end/2), real(GC_botsin.E_z) );
colormap('redbluehilight');
set( gca, 'ydir', 'normal' );
axis image;
% caxis( [-1,1] .* max(abs(real(GC_botsin.E_z(:)))) );
clim( 0.001 .* [-1,1] );
% caxis( 0.001 .* [-1,1] );
ylim( [-1.5, 1.5] );
% superimpose index contour
hold on;
contour( x, y - y(end/2), N, [2, 3.47], 'color', 'k', 'LineWidth', 1 );
colorbar;

% calc 1/e decay len
decay_len_botsin_mm = 1./imag(GC_botsin.k) * 1e-6











    
    
    

