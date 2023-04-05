% simulate aim cell for sopaipilla
% plots the decay length as a function of bottom duty cycle

clear; close all;

% dependencies
% laptop
addpath( genpath( [ 'C:\Users\howar\OneDrive\Desktop\Documents' ...
                    '\Boston University\Silicon Photonics\Grating Couplers' ...
                    '\grating_synthesis'] ) );
addpath( (genpath( [ 'C:\Users\howar\OneDrive\Desktop\Documents' ...
                     '\Boston University\Silicon Photonics\Grating Couplers' ...
                     '\grating_designs\aim_2021_07_16_SOPAs'])));

% simulate and plot for a bottom nitride grating

% set grating cell parameters
dxy             = 5;
lambda          = 1550;
k0              = 2*pi/lambda;
y_domain_size   = 4e3;

% set grating solver settings
num_modes   = 3;
BC          = 0;    % 0 = PEC
pml_options = [1, 100, 20, 2];
OPTS        = struct();

neff_slab   = 2.847003121289809;
n_clad      = 1.45;
theta       = 0; % deg
period      = 2*pi/( k0 * ( neff_slab - n_clad*sin( theta*pi/180 ) ) );
period      = round(period/dxy) * dxy;
guessk      = neff_slab*k0;

duty_cycles = 0.04 : dxy/period : 0.96;

% primary quantities to calculate
alpha   = zeros(size(duty_cycles));
beta    = zeros(size(duty_cycles));

% loop through all duty cycles
tic;
parfor duty_cycle = 1:length(duty_cycles)
    fprintf('(%d/%d) Running simulation for duty cycle = %g ...\n', duty_cycle, length(duty_cycles), duty_cycles(duty_cycle));
    
    % create grating cells
    GC = f_makeGratingCell_AIM( dxy, lambda, y_domain_size, period, 0, duty_cycles(duty_cycle), 0, OPTS);
    GC.numcells = 5;
    
    % run simulation
    GC = GC.runSimulation( num_modes, BC, pml_options, k0, guessk);

    % calculate alpha [nm^-1]
    alpha(duty_cycle) = imag(GC.k);

    % calculate beta [nm^-1]
    beta(duty_cycle) = real(GC.k);

end
toc;

% calculate 1/e^2 power decay length
decay_len_mm = 1 ./ alpha .* 1e-6 ./ 2;
neff = beta .* k0;

%%
close all;

min_dc = 100 * 150/545;
max_dc = 100 * (1 - (100/545));

% interpolate
duty_cycles_new     = 0.1:0.001:0.9;
alpha_new           = interp1(duty_cycles, alpha, duty_cycles_new);
decay_len_mm_new    = 1 ./ alpha_new .*1e-6 ./ 2;
beta_new            = interp1(duty_cycles, beta, duty_cycles_new);
neff_new            = beta_new ./ k0;

% add selected duty cycles for tapeout
x = [0.275, 0.45, 0.70, 0.75, 0.78, 0.798, 0.816];
y = zeros(size(x));
for idx = 1:length(x)
    idx2 = find(abs(x(idx) - duty_cycles_new) < 0.0001);
    y(idx) = decay_len_mm_new(idx2);
end

% % plot scattering strength (alpha) as a function of bottom duty cycle
% figure;
% plot(100*duty_cycles_new, alpha_new * 1e9, 'LineWidth', 1.5, 'Color', 'k');
% title( ['Scattering Strength vs. Bottom Duty Cycle @ ' ...
%         num2str(theta) '\circ and \lambda = ' ...
%         num2str(lambda) ' \mum'] );
% xlabel('Bottom Duty Cycle [%]'); ylabel('\alpha [m^{-1}]');
% %xlim([16 84]);
% grid on;

% plot field decay length as a function of bottom duty cycle
figure;
% plot(100*duty_cycles, decay_len_mm, 'LineWidth', 1.5, 'Color', 'r'); hold on;
plot(100*duty_cycles_new, decay_len_mm_new, 'LineWidth', 1.5, 'Color', 'b'); hold on;
plot([min_dc min_dc], [0 22], [max_dc max_dc], [0 22], 'LineWidth', 1.5, 'Color', 'k'); hold on;
scatter(100*x, y, 100, 'k', 'Marker', 'x', 'LineWidth', 2);
title( ['1/e^2 Power Decay Length vs. Bottom Duty Cycle @ ' ...
        num2str(theta) '\circ and \lambda = ' ...
        num2str(lambda) ' \mum'] );
xlabel('Bottom Duty Cycle [%]'); ylabel('Decay Length [mm]');
xlim([10 90]); ylim([0 20]);
grid on;

% plot effective index as a function of bottom duty cycle
% figure;
% plot(100*duty_cycles_new, neff, 'LineWidth', 1.5, 'Color', 'r');
% title( ['Effective Index vs. Bottom Duty Cycle @ ' ...
%         num2str(theta) '\circ and \lambda = ' ...
%         num2str(lambda) ' \mum'] );
% xlabel('Bottom Duty Cycle [%]'); ylabel('n_{eff}')
% grid on;









