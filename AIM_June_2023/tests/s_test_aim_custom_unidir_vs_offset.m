% custom aim gc cell simulation for SOPAIPILLA
% sweep through offset ratios to find maximum directionality

clear; close all;

% dependencies
addpath(genpath('C:\Users\howar\Siphot\Grating Couplers\grating_synthesis'));
addpath(genpath('C:\Users\howar\Siphot\Grating Couplers\grating_designs\AIM_June_2023'));

% simulation settings
num_modes       = 3;
BC              = 0;    % 0 = PEC
pml_options     = [1, 100, 20, 2];
OPTS            = struct();

% geometry parameters
thick_Si = 220;
space_Si_SiN = 100;
thick_bot_SiN = 60;
space_SiN_SiN = 0;
thick_top_SiN = 160;

geometry = struct('thick_Si', thick_Si, ...
                    'space_Si_SiN', space_Si_SiN, ...
                    'thick_bot_SiN', thick_bot_SiN, ...
                    'space_SiN_SiN', space_SiN_SiN, ...
                    'thick_top_SiN', thick_top_SiN ...
                    );

% duty cycles
dc_bot = 0.50;
dc_top = 0.50;

% cell parameters
dxy = 10;
lambda = 1550;
k0 = 2*pi/lambda;
y_domain_size = 4e3;

neff_slab = 2.842431701;
guess_k = neff_slab*k0;
period = 540;

% sweep offset ratios
offset_ratios = 0.4 : dxy/period : 0.42;

% variables to calculate
directivities = zeros(size(offset_ratios));
alpha = zeros(size(offset_ratios));

tic;
for i_off = 1:length(offset_ratios)
    fprintf('%d/%d Simulating for offset ratio = %g ...\n', i_off, length(offset_ratios), offset_ratios(i_off))

    % create grating cell and simulate
    GC = f_makeGratingCell_AIM_custom(dxy, lambda, y_domain_size, period, dc_top, dc_bot, offset_ratios(i_off), geometry);
    GC.numcells = 5;
    GC = GC.runSimulation(num_modes, BC, pml_options, k0, guess_k, OPTS);

    % calculate variables
    directivities(i_off) = 10*log10(GC.directivity);
    alpha(i_off) = imag(GC.k);

end
toc;

% 1/e^2 power decay length [mm]
decay_lens = 1 ./ alpha ./ 2 .* 1e-6;

% offsets [nm]
offsets = offset_ratios * period;

% plot results
figure;
plot(offsets, directivities, 'LineWidth', 1.5, 'Color', 'k');
title('Directivity vs. Offset Ratio');
xlabel('Offset Ratio');
ylabel('Directivity [dB]');
grid on;

figure;
plot(offsets, decay_lens, 'LineWidth', 1.5, 'Color', 'k');
title('1/e^2 Power Decay Length vs. Offset Ratio');
xlabel('Offset Ratio');
ylabel('Decay Length [mm]');
grid on;








