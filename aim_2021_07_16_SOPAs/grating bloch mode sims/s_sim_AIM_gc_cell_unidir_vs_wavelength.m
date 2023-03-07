% simulate aim cell for DII
% simulating unidirectionality as function of wavelength

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
                    '\grating_designs\aim_2021_07_16_SOPAs\bash_scripts' ...
                    '\2023_02_17_17_04_37_aimgc_geometry_b_lambda1550_optangle15_dx_5_up_clad_1d45'];
synth_obj       = load( [ synth_obj_path filesep 'synth_obj_aimgc.mat' ] );
synth_obj       = synth_obj.synth_obj;

% plotting decay length
% scatter_str_v_fills = synth_obj.sweep_variables.scatter_str_vs_fills;
% 
% decay_len_nm = 1./scatter_str_v_fills;
% decay_len_mm = decay_len_nm * 1e-6;
% 
% figure;
% imagesc( synth_obj.sweep_variables.fill_tops, synth_obj.sweep_variables.fill_bots, decay_len_mm );
% colorbar;
% xlabel('top duty cycle'); ylabel('bottom duty cycle');
% set(gca, 'ydir', 'normal');
% title('1/e field decay length (mm)');

% fills to use
fill_top = 0.34;
fill_bot = 0.46;

% pick out the rest of the parameters
indx_top        = find( fill_top == synth_obj.sweep_variables.fill_tops );
indx_bot        = find( fill_bot == synth_obj.sweep_variables.fill_bots );
period          = synth_obj.sweep_variables.periods_vs_fills( indx_bot, indx_top );
offset_ratio    = synth_obj.sweep_variables.offsets_vs_fills( indx_bot, indx_top ) ./ period;
k               = synth_obj.sweep_variables.k_vs_fills( indx_bot, indx_top );
%k               = 2*pi / synth_obj.lambda;
unidir          = synth_obj.sweep_variables.directivities_vs_fills( indx_bot, indx_top );
angle1550       = synth_obj.sweep_variables.angles_vs_fills( indx_bot, indx_top );

y_domain_size       = 4e3;
si_to_bot_sin       = 100;
bot_sin_thick       = 60;
sin_to_sin_thick    = 0;
top_sin_thick       = 160;

OPTS = struct( 'geometry','default si custom sin gratings full etch', ...
               'si_to_bot_sin', si_to_bot_sin, ...
               'bot_sin_thick', bot_sin_thick, ...
               'sin_to_sin_thick', sin_to_sin_thick, ...
               'top_sin_thick', top_sin_thick );

% set grating solver settings
num_modes   = 5;
BC          = 0;    % 0 = PEC
pml_options = [1, 100, 20, 2]; 
OPTS        = struct();
dxy         = 5;

si_to_bot_sin       = 100;
bot_sin_thick       = 60;
sin_to_sin_thick    = 0;
top_sin_thick       = 160;

OPTS = struct( 'geometry','default si custom sin gratings full etch', ...
               'si_to_bot_sin', si_to_bot_sin, ...
               'bot_sin_thick', bot_sin_thick, ...
               'sin_to_sin_thick', sin_to_sin_thick, ...
               'top_sin_thick', top_sin_thick );

% wavelengths to solve for
lambdas = 1450:10:1650;
lambdas_high = lambdas( lambdas > synth_obj.lambda );
lambdas_low  = lambdas( lambdas < synth_obj.lambda );

% saving variables
unidir_vs_lambda_high   = zeros(size(lambdas_high));
k_vs_lambda_high        = zeros(size(lambdas_high));
angle_vs_lambda_high    = zeros(size(lambdas_high));
unidir_vs_lambda_low    = zeros(size(lambdas_low));
k_vs_lambda_low         = zeros(size(lambdas_low));
angle_vs_lambda_low     = zeros(size(lambdas_high));

% first solve for longer wavelengths
% initial guessk 
guessk = k;

tic;
for ii = 1:length(lambdas_high)
    
    fprintf('simulating longer wavelengths, loop %i of %i\n', ii, length(lambdas_high));
    
    % update wavelength
    guessk = guessk * synth_obj.lambda./lambdas_high(ii);
    
    % make grating cell
    GC = f_makeGratingCell_AIM( dxy, lambdas_high(ii), y_domain_size, ...
                            period, fill_top, fill_bot, offset_ratio, OPTS  );
    GC.numcells = 5;

    % run simulation
    GC = GC.runSimulation( num_modes, BC, pml_options, 2*pi/lambdas_high(ii), guessk, OPTS );

    % get directionality + k
    unidir_vs_lambda_high(ii) = GC.directivity;
    k_vs_lambda_high(ii)      = GC.k;
    angle_vs_lambda_high(ii)  = GC.max_angle_up;
    
    % update guessk
    guessk = GC.k;
    
    toc;
    
end

% solve for shorter wavelengths
% initial guessk 
guessk = k;

tic;
for ii = 1:length(lambdas_low)
    
    fprintf('simulating shorter wavelengths, loop %i of %i\n', ii, length(lambdas_low));
    
    % update wavelength
    guessk = guessk * synth_obj.lambda./lambdas_low(end - (ii - 1));
    
    % make grating cell
    GC = f_makeGratingCell_AIM( dxy, lambdas_low(end - (ii - 1)), ...
                                y_domain_size, period, fill_top, fill_bot, ...
                                offset_ratio, OPTS  );
    GC.numcells = 5;

    % run simulation
    GC = GC.runSimulation( num_modes, BC, pml_options, 2*pi/lambdas_low(end - (ii - 1)), guessk, OPTS );

    % get directionality + k
    unidir_vs_lambda_low(end - (ii - 1)) = GC.directivity;
    k_vs_lambda_low(end - (ii - 1))      = GC.k;
    angle_vs_lambda_low(end - (ii - 1))  = GC.max_angle_up;
    
    % update guessk
    guessk = GC.k;
    
    toc;
    
end
    
% combine all variables
all_unidir_vs_lambda = [ unidir_vs_lambda_low, unidir,  unidir_vs_lambda_high ];
all_angles_vs_lambda = [ angle_vs_lambda_low, angle1550,  angle_vs_lambda_high ];
all_k_vs_lambda      = [ k_vs_lambda_low, k,  k_vs_lambda_high ];

figure;
yyaxis left;
plot( lambdas, 10*log10(all_unidir_vs_lambda), '-' );
xlabel('Wavelength (nm)', 'fontsize', 8); ylabel('Power up/Power down (dB)', 'fontsize', 8);
ylim([-60,60]);
makeFigureNice();
yyaxis right;
plot( lambdas, all_angles_vs_lambda, '-' );
ylabel('Angle (degrees)');
makeFigureNice();
f = gcf;
f.Position = [100   375   447   325];

figure;
plot( lambdas, 1e-7*1./(imag(all_k_vs_lambda)), '-' );
xlabel('Wavelength (nm)'); ylabel('1/e^2 power decay length (cm)');
makeFigureNice();



    
    
    

