% simulate period and effective index vs wavelength for
% vertical emission with bottom nitride gratings

clear;

% dependencies
% laptop
addpath( genpath( [ 'C:\Users\howar\OneDrive\Desktop\Documents' ...
                    '\Boston University\Silicon Photonics\Grating Couplers' ...
                    '\grating_synthesis'] ) );
addpath( [  'C:\Users\howar\OneDrive\Desktop\Documents' ...
            '\Boston University\Silicon Photonics\Grating Couplers' ...
            '\grating_designs\aim_2021_07_16_SOPAs'] );

y_domain_size   = 4e3;
num_modes       = 5;
BC              = 0;    % 0 = PEC
pml_options     = [1, 100, 20, 2];
OPTS            = struct();

% simulate and plot for a bottom nitride grating
duty_cycle = 0.75;

dxy     = 10;
lambda  = 1450:5:1650;
k0      = 2*pi ./ lambda;

n_clad      = 1.45;
theta       = 0; % deg

neff_init  = 2.84;
neff_new    = zeros(1, length(lambda));

period  = 2*pi ./ ( k0 .* ( neff_init - n_clad*sin( theta*pi/180 ) ) );
period_sim  = round(period ./ dxy) * dxy;

tic;
for i_lam = 1:length(lambda)
    % go through all wavelengths
    GC_botsin = f_makeGratingCell_AIM(  dxy, lambda(i_lam), y_domain_size, ...
                                        period_sim(i_lam), 0, duty_cycle, ...
                                        0, OPTS );
    GC_botsin.numcells = 5;

    neff_init   = 2.84;
    guessk      = neff_init * k0(i_lam);

    while neff_init - neff_new(i_lam) > 1e-3
        % test if neff has converged

        neff_init = neff_new(i_lam);

        % run simulation and get new neff
        GC_botsin = f_makeGratingCell_AIM(  dxy, lambda(i_lam), y_domain_size, ...
                                            period_sim(i_lam), 0, duty_cycle, ...
                                            0, OPTS );
        GC_botsin = GC_botsin.runSimulation(num_modes, BC, pml_options, ...
                                            k0(i_lam), guessk, OPTS );
        beta = real(GC_botsin.k);
        neff_new(i_lam) = beta / k0(i_lam);

        % set new period for next run
        period(i_lam)       = 2*pi/( k0(i_lam) * ( neff_new(i_lam) - n_clad*sin( theta*pi/180 ) ) );
        period_sim(i_lam)   = round( period_sim(i_lam) / dxy )  * dxy;
        guessk              = neff_new(i_lam) * k0(i_lam);

    end
    toc;
end

theta = asin(neff_new - lambda ./ period) * 180/pi;

%%
close all;

figure;
yyaxis left;
plot(lambda, period, 'LineWidth', 1.5);
title('Period/n_{eff} vs. Wavelength for Vertical Emission w/ Bottom Nitride Gratings');
xlabel('Wavelength [nm]');
ylabel('Period [nm]');

yyaxis right;
plot(lambda, neff_new, 'LineWidth', 1.5);
ylabel('n_{eff}');

grid on;

figure;
plot(period, theta, 'LineWidth', 1.5, 'Color', 'k');
title(['Radiation Angle vs. Period @ ' num2str(duty_cycle*100) '%']);
xlabel('Period [nm]');
ylabel('\theta [\circ]');






