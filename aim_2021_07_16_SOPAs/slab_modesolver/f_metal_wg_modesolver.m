function [ field, k, x ] = f_metal_wg_modesolver( lambda0, n, layer_t, dx, guessk, n_modes )
% written by bohan zhang
% 
% 1D slab FDFD modesolver
% all units in um
% 
% Inputs:
%     lambda0
%       type: scalar, double
%       desc: wavelength in um
%     n
%       type: scalar, double
%       desc: refractive index of waveguide core
%     layer_t
%       type: scalar, double
%       desc: thickness of waveguide, in um
%     dx
%       type: scalar, double
%       desc: discretization in um
%     guessk
%       type: scalar, double
%       desc: guess propagation value, units rad/um
%     n_modes
%       type: scalar, integer
%       desc: number of modes to (try to) solve for
% 
% Outputs
%     field
%       type: matrix, double
%       desc: field, dimensions of field(x) vs. mode #
%     k
%       type: vector, double
%       desc: propagation constant in y (direction of prop) vs. mode #
%     x
%       type: vector, double
%       desc: x coordinate vector


% constants
ep_0        = 8.854e-12 * 1e-6;             % units F/um
mu_0        = 4*pi*1e-07 * 1e-6;            % units H/umm
c           = 3e8 * 1e6;                    % units um/s
omega       = 2*pi*c/lambda0;               % units rad/s
ep_r        = n^2;                          % unitless


% space coordinate vector
x = 0:dx:layer_t;

% check space is discretized well
if abs( x(end) - layer_t ) > 1e-6
    fprintf('Warning: space is not integer discretized\n'); 
end


% generate dx^2 operator
nx      = length(x);
dx2_op  = (1/(dx^2)) * ( diag( -2*ones(nx, 1), 0 ) + diag( ones(nx-1, 1), 1 ) + diag( ones(nx-1, 1), -1 ) );
dx2_op( 1, 1 )      = 0;
dx2_op( end, end )  = 0;

% guess eigenvalue
guesskx_sq = omega^2*mu_0*ep_0*ep_r - guessk^2;

% solve eigs
[field, kx_sq_diag] = eigs(dx2_op, n_modes, guesskx_sq);
% flip signs of stuff
kx_sq_diag          = -kx_sq_diag;
field               = -field;

% calculate propagation constants
k  = sqrt( omega^2*mu_0*ep_0*ep_r - diag(kx_sq_diag) );

end



