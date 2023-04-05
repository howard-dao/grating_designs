% solves TE mode of waveguide with PEC boundaries

clear; close all;

% --------------------------------------------------------------------|
% Parameter def.
% --------------------------------------------------------------------|

% wave parameters
lambda0     = 1.55;                                 % units um
k0          = 2*pi/lambda0;

% define material layers
n_core = 1.5;

% discretization
dx = 0.001;                                          % units um

% TEMP calculating cutoff d
% constants
ep_0        = 8.854e-12 * 1e-6;             % units F/um
mu_0        = 4*pi*1e-07 * 1e-6;            % units H/umm
c           = 3e8 * 1e6;                    % units um/s
omega       = 2*pi*c/lambda0;               % units rad/s
ep_r        = n_core^2;                     % unitless
d_cutoff    = 2*pi/( omega * sqrt( mu_0 * ep_0 * ep_r )  );


% define layer thicknesses
core_d = 2.0;                                   % core thickness, units um

% calculate analytical k of the fundamental mode
k_fundamental = sqrt( omega^2*mu_0*ep_0*ep_r - 2*pi/core_d ); 

% guess k
guessk = k_fundamental;                             % units rad/um

% choose # of modes to solve for
n_modes = 2;

% use modesolver function
field_all = struct( 'field', {}, 'k', {} );  % saves all the fields

for ii = 1:length(dx) % trying different discretizations
    
    [ field, k, x ] = f_metal_wg_modesolver( lambda0, n_core, core_d, dx, guessk, n_modes );
    
    field_all(ii).field = field;
    field_all(ii).k     = k;
    
end

% % compute analytical solution (symmetric only)
% [ neff, k_analytical, kx_temp, alpha_temp, field_temp ] = solve_symm_slab( core_d*1e-6, n_core, n_core, lambda0*1e-6, 'TE', true );

% % plot the index
% figure;
% plot(x, n_array);
% title('Index distribution');
% xlabel('x (m)'); ylabel('Index');
% makeFigureNice();

% plot the modes
for i_mode = 1:size(field,2)

        % plots index and field
%     figure;
%     plotyy( x, (field(:,i_mode)), x, n_array ); hold on;
%     makeFigureNice();
%     title(sprintf('mode: %i of %i', i_mode, size(field,2)));
%     xlabel('x');
%     legend('Field', 'index');
%     makeFigureNice();

    % plot field only
    figure;
    plot( x, (field(:,i_mode)) ); hold on;
    title(sprintf('mode: %i of %i', i_mode, size(field,2)));
    xlabel('x');  ylabel('E field');
    makeFigureNice();
    
end

% 
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








