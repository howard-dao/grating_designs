% simulate aim cell for sopaipilla
% plots directionality and 1/e^2 decay length as a function of offset ratio

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

% set grating cell parameters
dxy = 5;
y_domain_size   = 4e3;

% fills to use
fill_top = 0.50;
fill_bot = 0.74;

% pick out the rest of the parameters
indx_top        = find( abs(fill_top - synth_obj.sweep_variables.fill_tops) < 0.001 );
indx_bot        = find( abs(fill_bot - synth_obj.sweep_variables.fill_bots) < 0.001 );
% period          = synth_obj.sweep_variables.periods_vs_fills( indx_bot, indx_top );
period = 550;
% offsets         = 0 : synth_obj.discretization : period;
offsets         = 0 : dxy : period;
offset_ratio    = offsets ./ period;
k               = synth_obj.sweep_variables.k_vs_fills( indx_bot, indx_top );

% set grating solver settings
num_modes   = 2;
BC          = 0;    % 0 = PEC
pml_options = [1, 100, 20, 2];
OPTS        = struct();
guessk      = k;

directivity = zeros(size(offset_ratio));
alpha       = zeros(size(offset_ratio));

tic;
parfor i_offset = 1:length(offset_ratio)
    fprintf('%d/%d Running simulation for lower duty cycle = %g, upper duty cycle = %g, and offset ratio = %g ...\n', ...
            i_offset, length(offset_ratio), fill_top, fill_bot, offset_ratio(i_offset));

    % create grating cells
    GC = f_makeGratingCell_AIM( dxy, synth_obj.lambda, y_domain_size, ...
                                period, fill_top, fill_bot, ...
                                offset_ratio(i_offset), OPTS  );
    GC.numcells = 5;
    
    % run simulation
    GC = GC.runSimulation( num_modes, BC, pml_options, 2*pi/synth_obj.lambda, guessk);
    
    % calculate directionality [dB]
    directivity(i_offset) = 10*log10(GC.directivity);

    % calculate alpha [nm^-1]
    alpha(i_offset) = imag(GC.k);

end
toc;

% 1/e^2 power decay length [mm]
decay_len_mm = 1 ./ alpha .* 1e-6 ./ 2;

%%
close all;

figure;
plot(offset_ratio, directivity, 'LineWidth', 1.5, 'Color', 'k');
title(['Directionality vs. Offset Ratio: D_l=' num2str(fill_bot*100) '%, D_u=' num2str(fill_top*100) '%']);
xlabel('Offset Ratio');
ylabel('Directionality [dB]');
grid on;

figure;
plot(offset_ratio, decay_len_mm, 'LineWidth', 1.5, 'Color', 'k');
title(['Decay Length vs. Offset Ratio: D_l=' num2str(fill_bot*100) '%, D_u=' num2str(fill_top*100) '%']);
xlabel('Offset Ratio');
ylabel('Decay Length [mm]');
ylim([0 inf]);
grid on;











    
    
    

