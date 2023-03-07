% written by Bohan Zhang
% 
% Simulating first mode of waveguide, comparing for different wavelengths

clc; clear; close all;

% --------------------------------------------------------------------|
% Parameter def.
% --------------------------------------------------------------------|

% wave parameters
lambda0     = 1.50e-6;    % m
k0          = 2*pi/lambda0;
dlambda     = 0.1e-6;
lambda      = (lambda0 - dlambda):dlambda:(lambda0 + dlambda); % all wavelengths to simulate
lambda      = lambda0;
k           = 2*pi./lambda;

% define material layers
n_clad      = 1.0;
n_core      = 3.5;
n           = [ n_clad, n_core, n_clad ];

% discretization
dx  = 0.1*lambda0;

% define layer thicknesses
clad_d    = 5*lambda0;  % layer depths
core_d    = 1*lambda0; % layer depths
layer_t   = [ clad_d, core_d, clad_d ];

% choose # of modes to solve for
n_modes     = 1;

% use modesolver function
field_all = struct( 'field', {}, 'k', {} );  % saves all the fields

% simulate vs wavelength
for ii = 1:length(k)
    
    guessk = k(ii);
    
    [ field, k_fdfd, x, n_array ]  = slab_modesolver_bz( lambda(ii), n, layer_t, dx, guessk, n_modes );
    
    field_all(ii).field = field;
    field_all(ii).k     = k_fdfd;
    
end

% compute analytical solution (symmetric only)
[ neff, k_analytical, kx_temp, alpha_temp, field_temp ] = solve_symm_slab( core_d, n_clad, n_core, lambda0, 'TE', false );

% plot the index
figure;
plot(x, n_array);
title('Index distribution');
xlabel('x (m)'); ylabel('Index');
makeFigureNice();

% plot the modes vs wavelength
figure;
for ii = 1:length(field_all)

    cur_field = field_all(ii).field;
   
    plotyy( x, cur_field, x, n_array ); hold on;
    makeFigureNice();
%     title(sprintf('mode: %i of %i', ii, size(field,2)));
    xlabel('x');
%     legend('Field', 'index');
    makeFigureNice();
    
end

% % plot analytical beta versus FDFD beta
% figure;
% plot( k_analytical, '--o' ); hold on;
% legendstrs = {};
% for ii = 1:length(field_all)
%     plot( field_all(ii).k, '--x' );
%     legendstrs{end+1} = sprintf( 'FDFD, dx = %f lambda0', dx(ii)/lambda0 );
% end
% ylabel(' \beta aka k_z');
% legend( { 'Analytical', legendstrs{:} } );
% title('Analytical vs FDFD solved prop. constants');
% makeFigureNice();








