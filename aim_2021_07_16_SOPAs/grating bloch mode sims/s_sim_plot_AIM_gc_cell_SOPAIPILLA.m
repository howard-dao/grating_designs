% simulate aim cell for DII

clear; close all;

% dependencies
% laptop
addpath( genpath( [ 'C:\Users\howar\OneDrive\Desktop\Documents' ...
                    '\Boston University\Silicon Photonics\Grating Couplers' ...
                    '\grating_synthesis'] ) );
addpath( [  'C:\Users\howar\OneDrive\Desktop\Documents' ...
            '\Boston University\Silicon Photonics\Grating Couplers' ...
            '\grating_designs\aim_2021_07_16_SOPAs'] );

% load unidir data
synth_obj_path  = [ 'C:\Users\howar\OneDrive\Desktop\Documents' ...
                    '\Boston University\Silicon Photonics\Grating Couplers' ...
                    '\nitride_space100_top_nitride220_lam1550_angle0_dx_5_up'];
synth_obj       = load( [ synth_obj_path filesep 'synth_obj_aimgc.mat' ] );
synth_obj       = synth_obj.synth_obj;

% fills to use
fill_top = 0.50;
fill_bot = 0.74;

% pick out the rest of the parameters
indx_top        = find( abs(fill_top - synth_obj.sweep_variables.fill_tops) < 0.001 );
indx_bot        = find( abs(fill_bot - synth_obj.sweep_variables.fill_bots) < 0.001 );
period          = synth_obj.sweep_variables.periods_vs_fills( indx_bot, indx_top );
% offset          = synth_obj.sweep_variables.offsets_vs_fills( indx_bot, indx_top );
% offset_ratio    = offset ./ period;
% offset_ratio = 0.376147;
offset_ratio = 0.165;
k               = synth_obj.sweep_variables.k_vs_fills( indx_bot, indx_top );

y_domain_size       = 4e3;
si_to_bot_sin       = 100;
bot_sin_thick       = 220;
sin_to_sin_thick    = 100;
top_sin_thick       = 220;

OPTS = struct( 'geometry',          'default si custom sin gratings full etch', ...
               'si_to_bot_sin',     si_to_bot_sin, ...
               'bot_sin_thick',     bot_sin_thick, ...
               'sin_to_sin_thick',  sin_to_sin_thick, ...
               'top_sin_thick',     top_sin_thick );

GC = f_makeGratingCell_AIM( 5, synth_obj.lambda, y_domain_size, ...
                            period, fill_top, fill_bot, ...
                            offset_ratio, OPTS  );
GC.numcells = 5;

% plot index
% GC.plotIndex();

% set grating solver settings
num_modes   = 5;
BC          = 0;    % 0 = PEC
pml_options = [1, 100, 20, 2];
guessk      = k;
OPTS        = struct();

% run simulation
tic;
GC = GC.runSimulation( num_modes, BC, pml_options, 2*pi/synth_obj.lambda, guessk);%, OPTS );
toc;

% GC.directivity = 10*log10(GC.directivity);

% plot e field
% GC.plot_E_field_gui();

% plot field with waveguide overlay
x = 1e-3 * ( 0 : GC.dx : GC.domain_size(2)*GC.numcells - GC.dx );
y = GC.y_coords*1e-3;
N = repmat( GC.N, 1, GC.numcells );

figure;
imagesc( x, y - y(end/2), real(GC.E_z) );
colormap('redbluehilight');
set( gca, 'ydir', 'normal' );
% set(gca, 'XTick', [], 'YTick', []);
axis image;
clim( 0.001 .* [-1,1] );
ylim( [-1.5, 1.5] );
% superimpose index contour
hold on;
contour( x, y - y(end/2), N, [2, 3.47], 'color', 'k', 'LineWidth', 1 );
colorbar;

% calc 1/e decay len
decay_len_mm = 1./imag(GC.k) * 1e-6
    
%% Simulate single layer gratings
% simulate and plot for a bottom nitride grating of similar duty cycle

dxy     = 4;
lambda  = 1550;
k0      = 2*pi/lambda;

% neff_slab   = 2.842431701;
neff_slab   = 2.847003121289809;
n_clad      = 1.45;
theta       = 0; % deg
period      = 2*pi/( k0 * ( neff_slab - n_clad*sin( theta*pi/180 ) ) );
period      = round(period/dxy)  * dxy;
guessk      = neff_slab*k0;

fill_bot = 0.90;

GC_botsin = f_makeGratingCell_AIM( dxy, lambda, y_domain_size, period, 0, fill_bot, 0, OPTS );
GC_botsin.numcells = 5;
                           
% plot index
% GC_botsin.plotIndex();

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
set(gca, 'XTick', [], 'YTick', []);
axis image;
clim( 0.0005 .* [-1,1] );
ylim( [-1.5, 1.5] );
% superimpose index contour
hold on;
contour( x, y - y(end/2), N, [2, 3.47], 'color', 'k', 'LineWidth', 1 );
colorbar;

% calc 1/e decay len
decay_len_botsin_mm = 1./imag(GC_botsin.k) * 1e-6




% simulate and plot for a upper nitride grating of similar duty cycle
% duty_cycle = 0.25;

dxy     = 10;
lambda  = 1550;
k0      = 2*pi/lambda;

neff_slab   = 2.842431701;
n_clad      = 1.45;
% theta       = 0; % deg
period      = 2*pi/( k0 * ( neff_slab - n_clad*sin( theta*pi/180 ) ) );
period      = round(period/dxy)  * dxy;
guessk      = neff_slab*k0;

GC_topsin = f_makeGratingCell_AIM(  dxy, lambda, y_domain_size, ...
                                    period, fill_top, 0, 0, OPTS );
GC_topsin.numcells = 5;
                           
% plot index
% GC_topsin.plotIndex();

% run simulation
tic;
GC_topsin = GC_topsin.runSimulation( num_modes, BC, pml_options, k0, guessk, OPTS );
toc;

% plot e field
% GC.plot_E_field_gui();

% plot field with waveguide overlay
x = 1e-3 * ( 0 : GC_topsin.dx : GC_topsin.domain_size(2)*GC_topsin.numcells - GC_topsin.dx );
y = GC_topsin.y_coords*1e-3;
N = repmat( GC_topsin.N, 1, GC_topsin.numcells );

figure;
imagesc( x, y - y(end/2), real(GC_topsin.E_z) );
colormap('redbluehilight');
set( gca, 'ydir', 'normal' );
axis image;
clim( 0.001 .* [-1,1] );
ylim( [-1.5, 1.5] );
% superimpose index contour
hold on;
contour( x, y - y(end/2), N, [2, 3.47], 'color', 'k', 'LineWidth', 1 );
colorbar;

% calc 1/e decay len
decay_len_topsin_mm = 1./imag(GC_topsin.k) * 1e-6













    
    
    

