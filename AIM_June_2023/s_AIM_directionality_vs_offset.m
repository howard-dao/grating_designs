% testing AIM grating cell

clear; close all;

% dependencies
addpath(genpath([   'C:\Users\howar\OneDrive\Desktop\Documents\Boston University' ...
                    '\Silicon Photonics\Grating Couplers\grating_synthesis']));

dxy     = 10;
lambda  = 1550;
k0      = 2*pi/lambda;

neff_slab   = 2.842431701;
n_clad      = 1.45;
theta       = 15; % deg
period      = 2*pi/( k0 * ( neff_slab - n_clad*sin( theta*pi/180 ) ) );
period      = round(period/dxy)  * dxy;
y_domain_size = 4e3;

fill_top = 0.9;
fill_bot = 0.9;

offset_ratios = linspace(0,0.99, 10);

% saving variables
directionality_vs_offset = zeros(size(offset_ratios));

% set grating solver settings
num_modes   = 5;
BC          = 0;    % 0 = PEC
pml_options = [1, 100, 20, 2];
OPTS        = struct();
guessk      = k0 * neff_slab;

tic;
for i_offset = 1:length(offset_ratios)
    
    fprintf('loop %i of %i\n', i_offset, length(offset_ratios));

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

    GC = f_makeGratingCell_AIM_custom(  dxy, lambda, y_domain_size, period, fill_top, fill_bot, ...
                                        offset_ratios(i_offset), geometry );

    % run simulation
    GC = GC.runSimulation( num_modes, BC, pml_options, k0, guessk, OPTS );
    
    % update k
    guessk = GC.k;
    
    % save directionality
    directionality_vs_offset(i_offset) = GC.directivity;
    
    toc;
    
end

% plot e field
% GC.plot_E_field_gui();

% directionality vs. offset
figure;
plot( offset_ratios, directionality_vs_offset, '-o' );
xlabel('offset ratio'); ylabel('directionality up/down');
makeFigureNice();

















