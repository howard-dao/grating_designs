% written by Bohan Zhang
% 
% created: 11/4/16 6:43 pm
% last updated: 1/26/17
% 
% testing Slab Waveguide Solver, using Finite Difference technique

clear; close all;

% --------------------------------------------------------------------|
% Parameter def.
% --------------------------------------------------------------------|

% wave parameters
lambda0     = 1.55e-6;    % m
k0          = 2*pi/lambda0;

% define material layers
n_clad      = 1.44;
% n_core      = 3.477;
n_core      = 1.44;
n           = [ n_clad, n_core, n_clad ];

% discretization
% dx          = [ 0.1*lambda0, 0.05*lambda0, 0.01*lambda0 ];
dx          = 0.05*1e-6;

% define layer thicknesses
clad_d    = 5*1e-6;  % layer depths
core_d    = 0.15*1e-6; % layer depths
layer_t     = [ clad_d, core_d, clad_d ];

% guess k
guessk      = k0*n_core;

% choose # of modes to solve for
n_modes     = 1;

% use modesolver function
field_all = struct( 'field', {}, 'k', {} );  % saves all the fields

for ii = 1:length(dx) % trying different discretizations
    
    [ field, k_fdfd, x, n_array ]  = slab_modesolver_bz( lambda0, n, layer_t, dx(ii), guessk, n_modes );
    
    field_all(ii).field = field;
    field_all(ii).k     = k_fdfd;
    
end

% compute analytical solution (symmetric only)
[ neff, k_analytical, kx_temp, alpha_temp, field_temp ] = solve_symm_slab( core_d, n_clad, n_core, lambda0, 'TE', true );

% % plot the index
% figure;
% plot(x, n_array);
% title('Index distribution');
% xlabel('x (m)'); ylabel('Index');
% makeFigureNice();

% plot the modes
for i_mode = 1:size(field,2)

    figure;
    plotyy( x, (field(:,i_mode)), x, n_array ); hold on;
    makeFigureNice();
    title(sprintf('mode: %i of %i', i_mode, size(field,2)));
    xlabel('x');
    legend('Field', 'index');
    makeFigureNice();
    
end


% plot analytical beta versus FDFD beta
figure;
plot( k_analytical, '--o' ); hold on;
legendstrs = {};
for ii = 1:length(field_all)
    plot( field_all(ii).k, '--x' );
    legendstrs{end+1} = sprintf( 'FDFD, dx = %f lambda0', dx(ii)/lambda0 );
end
ylabel(' \beta aka k_z');
legend( { 'Analytical', legendstrs{:} } );
title('Analytical vs FDFD solved prop. constants');
makeFigureNice();








