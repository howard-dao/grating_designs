% simulate aim cell for DII

clear;

% dependencies
% laptop
addpath( genpath( [ 'C:\Users\howar\OneDrive\Desktop\Documents' ...
                    '\Boston University\Silicon Photonics\Grating Couplers' ...
                    '\grating_synthesis'] ) );
addpath( (genpath( [ 'C:\Users\howar\OneDrive\Desktop\Documents' ...
                     '\Boston University\Silicon Photonics\Grating Couplers' ...
                     '\grating_designs\aim_2021_07_16_SOPAs'])));

% simulate and plot for a bottom nitride grating of similar duty cycle
duty_cycles = 0.1:0.01:0.9;

% set grating cell parameters
dxy             = 10;
lambda          = 1550;
y_domain_size   = 4e3;
k0              = 2*pi/lambda;
OPTS            = struct();

% set grating solver settings
num_modes   = 5;
BC          = 0;    % 0 = PEC
pml_options = [1, 100, 20, 2];

neff_slab   = 2.842431701;
n_clad      = 1.45;
theta       = 0; % deg
period      = 2*pi/( k0 * ( neff_slab - n_clad*sin( theta*pi/180 ) ) );
period      = round(period/dxy)  * dxy;
guessk      = neff_slab*k0;

alpha           = zeros(length(duty_cycles), 1);
decay_len_mm    = zeros(length(duty_cycles), 1);
beta            = zeros(length(duty_cycles), 1);
neff            = zeros(length(duty_cycles), 1);

% loop through all duty cycles
tic;
for ii = 1:length(duty_cycles)
    GC = f_makeGratingCell_AIM( dxy, lambda, y_domain_size, ...
                                period, 0, duty_cycles(ii), 0, OPTS );
    GC.numcells = 5;
    
    % run simulation
    GC = GC.runSimulation( num_modes, BC, pml_options, k0, guessk, OPTS );

    % calculate alpha and 1/e decay length
    alpha(ii)           = imag(GC.k);
    decay_len_mm(ii)    = 1 / alpha(ii) * 1e-6;
    
    % calculate beta and effective index
    beta(ii)            = real(GC.k);
    neff(ii)            = beta(ii) / k0;
end
toc;
%%
close all;

% plot scattering strength (alpha) as a function of bottom duty cycle
figure;
plot(duty_cycles, alpha * 1e9, 'LineWidth', 1.5, 'Color', 'k');
title( ['Scattering Strength vs. Bottom Duty Cycle @ ' ...
        num2str(theta) '\circ and \lambda = ' ...
        num2str(lambda) ' \mum'] );
xlabel('Bottom Duty Cycle'); ylabel('\alpha [m^{-1}]');
%xlim([0.16 0.84]);
grid on;

% plot decay length as a function of bottom duty cycle
figure;
plot(duty_cycles, decay_len_mm, 'LineWidth', 1.5, 'Color', 'b');
title( ['Decay Length vs. Bottom Duty Cycle @ ' ...
        num2str(theta) '\circ and \lambda = ' ...
        num2str(lambda) ' \mum'] );
xlabel('Bottom Duty Cycle'); ylabel('Decay Length [mm]');
xlim([0.16 0.84]); ylim([0 inf]);
grid on;

% plot effective index as a function of bottom duty cycle
figure;
plot(duty_cycles, neff, 'LineWidth', 1.5, 'Color', 'r');
title( ['Effective Index vs. Bottom Duty Cycle @ ' ...
        num2str(theta) '\circ and \lambda = ' ...
        num2str(lambda) ' \mum'] );
xlabel('Bottom Duty Cycle'); ylabel('n_{eff}')
grid on;


