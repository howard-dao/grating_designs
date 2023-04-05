function [ field, beta, x, n_array ] = f_slab_modesolver_bz( lambda0, n, layer_t, dx, guessk, n_modes, PML_options )
% written by bohan zhang
% 
% 1D slab FDFD modesolver
% all units in um
%
% actually i believe method is unit-less but whatever
%
% also I think this only solves for TE 
% 
% Inputs:
%     lambda0
%         : wavelength (um)
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
%     PML_options:
%       type: struct
%       desc: PML options
%               PML_options(1): PML in x direction (yes=1 or no=0)
%               PML_options(2): length of PML layer, unit 
%               PML_options(3): strength of PML in the complex plane
%               PML_options(4): PML polynomial order (1, 2, 3...)              
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
ep_0        = 8.854e-12 * 1e-6;             % units F/um
mu_0        = 4*pi*1e-07 * 1e-6;            % units H/um
c           = 3e8 * 1e6;                    % units um/s
omega       = 2*pi*c/lambda0;               % units rad/s

% calculate number of discrete points
Nx      = layer_t/dx;
N_total = sum(Nx(:));
if is_integer_matrix( Nx ) == false         % testing if space is discretized well
   fprintf('Warning: space is not integer discretized\n'); 
end
Nx      = int16(Nx);
N_total = int16(N_total);
% Nx      = floor(Nx);

% space coordinate vectorj
x = linspace(0, sum(layer_t(:)), N_total);

% make material
ep_r            = n.^2;
ep_r_array      = [];                                       % index array
for ii = 1:length(ep_r)
    ep_r_array = [ ep_r_array, ep_r(ii)*ones(1,Nx(ii)) ];
end
ep_r_diag       = diag( ep_r_array );                       % convert to diagonal matrix
n_array         = sqrt( ep_r_array );

% make pmls
% draw in PMLs
if PML_options(1) == 1
    
    % grab params
    pml_len_um  = PML_options(2);   % length of pml in um
    pml_str     = PML_options(3);   % strength of pml in complex plane
    pml_order   = PML_options(4);   % pml polynomial order
    
%     ny_pml = 2 * round(pml_len_nm/disc);                                             % number of discretizations that pml spans, double sampled grid
%     if abs(ny_pml - round(ny_pml)) >= 1e-5
%         % discretization was not integer value
%         error('Integer # of discretizations did not fit into the PML');
%     end
%     y_indx = 1:N_total;
    
    % local coordinates for pml, on double grid
    nx_pml = 2*round(pml_len_um/dx);
    x_pml  = 1:nx_pml;
    
    % pml materials only
    pml_x   = (1 + 1i * pml_str * ( x_pml./nx_pml ).^( pml_order )).';
    
    % draw stretched coordinate pml
    pml_x_all                        = ones( 2*N_total, 1 );
    pml_x_all( 1:nx_pml )          = flipud( pml_x(1:end) );
    pml_x_all( end-nx_pml+1:end )    = pml_x;
    
    % stretched coordinate operator
%     pml_y_all_vec   = pml_y_all(:);
    Sx_f            = diag( 1./pml_x_all(2:2:end), 0);                % half step for forward Sy
    Sx_b            = diag( 1./pml_x_all(1:2:end-1), 0);              % on grid for backwards Sy
    
end

% make derivative operators
Dx_f            = (1/(dx))*( diag( -1*ones( size(ep_r_diag,1), 1), 0) + diag( 1*ones(size(ep_r_diag,1)-1,1), 1) );  % forward
Dx_b            = (1/(dx))*( diag( -1*ones( size(ep_r_diag,1)-1, 1), -1) + diag( 1*ones(size(ep_r_diag,1),1), 0) ); % backwards
% Dx_b(1,1)       = 2*Dx_b(1,1); % i think this is necessary to account for half grid
% Dx_f(1,2)       = 1/2*Dx_f(1,2);
% Dx_f(1,1)       = 2*Dx_f(1,1);
% Dx_b(end,end)       = 0; % i think this is necessary to account for half grid
% Dx_f(end,end)   = 2*Dx_f(end,end); % i think this is necessary to account for half grid

Dx2             = (1/(dx^2)) * ( diag( 1*ones( size(ep_r_diag,1)-1, 1 ), -1 ) ...
                + diag( 1*ones( size(ep_r_diag,1)-1, 1 ), 1 ) ...
                + diag( -2*ones( size(ep_r_diag,1), 1 ), 0 ) );



% apply pml operators
if PML_options(1) == 1
%     Dx_f2 = Sx_f * Dx_f;
%     Dx_b2 = Sx_b * Dx_b;
%     Dx2 = Dx_b2 * Dx_f2;
%     
%     Dx2_v6 = Sx_b * Sx_f * Dx_b * Dx_f;
%     
%     Dx2_v2 = Sx_b * Sx_f * Dx2_v1;
%     Dx2_v3 = Sx_f * Sx_b * Dx2_v1;
%     Dx2_v4 = Sx_b * Sx_b * Dx2_v1;
%     Dx2_v5 = Sx_f * Sx_f * Dx2_v1;
    
%     Dx2 = Sx_f * Sx_b * Dx2;
    
    % what if i try this
    ep_r_diag = diag( ep_r_array.' .* (1./pml_x_all(1:2:end-1)), 0 );
    
end



% Dx2(1,1)        = 0;
% Dx2(end,end)    = 0;
% Dx2(2,1)        = 0;
% Dx2(end-1,end)    = 0;
    

% % DEBUG plot dx f dx b
% figure;
% spy( Dx_f );
% % set(gca, 'ydir', 'normal');
% title('Dx_f');
% % colorbar;

% make eigen matrix to solve
% A =  Dx_b*Dx_f + (omega^2)*ep_0*(1/mu_0)*ep_r_diag;
A = (Dx2 + (omega^2)*ep_0*mu_0*ep_r_diag);

% solve eigs
[field, beta_sq_diag] = eigs( A, n_modes, guessk^2 );

% calculate propagation constants
beta = sqrt( diag(beta_sq_diag) );

% flip field signage for some reason
field = -field;

end

