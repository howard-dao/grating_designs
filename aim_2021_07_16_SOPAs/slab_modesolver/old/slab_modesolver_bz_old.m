function [ field, beta, x, n_array ] = slab_modesolver_bz( lambda0, n, layer_t, dx, guessk, n_modes )
% written by bohan zhang
% 
% 1D slab FDFD modesolver
% all units in meters (ideally units would be in um but whatever)
% 
% Inputs:
%     lambda0
%         : wavelength (m)
%     n
%         : array of indexes
%     layer_t
%         : thickness of each index layer
%     dx
%         : discretization
%     guessk
%         : guess value of propagation constant (2*pi*neff/lambda0)
%     n_modes
%         : number of modes to (try to) solve for
% 
% Outputs
%     field
%         : has dimensions of ( field(x), mode# )
%     beta
%         : propagation constant in z per mode
%     x
%         : x coordinate vector
%     n_array
%         : array of indices


% constants
ep_0        = 8.854e-12;
mu_0        = 4*pi*1e-07;   % H per m
c           = 3e8;        % m/s
omega       = 2*pi*c/lambda0;

% calculate number of discrete points
Nx      = layer_t/dx;
N_total = sum(Nx(:));
if is_integer_matrix( Nx ) == false         % testing if space is discretized well
   fprintf('Warning: space is not integer discretized\n'); 
end
Nx      = int16(Nx);
N_total = int16(N_total);
% Nx      = floor(Nx);

% space coordinate vector
x = linspace(0, sum(layer_t(:)), N_total);

% make material
ep_r            = n.^2;
ep_r_array      = []; % index array
for ii = 1:length(ep_r)
    ep_r_array = [ ep_r_array, ep_r(ii)*ones(1,Nx(ii)) ];
end
ep_r_diag       = diag( ep_r_array );      % convert to diagonal matrix
n_array         = sqrt( ep_r_array );

% make derivative operators
Dx_f            = (1/(dx))*( diag( -1*ones( size(ep_r_diag,1), 1), 0) + diag( 1*ones(size(ep_r_diag,1)-1,1), 1) );  % forward
Dx_b            = (1/(dx))*( diag( -1*ones( size(ep_r_diag,1)-1, 1), -1) + diag( 1*ones(size(ep_r_diag,1),1), 0) ); % backwards

% make eigen matrix to solve
A =  mu_0*(Dx_b*(1/mu_0)*Dx_f + (omega^2)*ep_0*ep_r_diag);

% solve eigs
[field, beta_sq_diag] = eigs(A, n_modes, guessk^2);

% calculate propagation constants
beta = sqrt( diag(beta_sq_diag) );

end

