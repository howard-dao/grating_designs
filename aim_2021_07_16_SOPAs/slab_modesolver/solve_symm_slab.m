function [ neff, beta, kx, alpha, Field ] = solve_symm_slab( d, ncl, nco, lambda0, pol, show_plot )
% Written by: Bohan Zhang (bhnzhang@gmail.com)
% 
% Last updated: 11/18/16 9:58am
% 
% Solves ANALYTICAL solution for symmetric slab waveguide
%     Based on Yariv - ch 3
% 
% solving for TE/TM modes
% 
% Inputs:
%     d
%         :   width of slab
%             slab spans -d/2 to d/2
%     ncl
%         :   index of cladding
%     nco
%         :   index of core
%     lambda0
%         :   wavelength
%     pol
%         :   string, 'TE' for TE, 'TM' for TM
%     show_plot
%         :   true/false, show effective index curve flag (OPTIONAL)
%             defaults to false
% 
% Outputs:
%     neff
%         :   effective index for each mode
%     beta
%         :   propagation constant for each mode
%       kx
%         :   transverse guided k for each mode
%       alpha
%         :   transverse exponential decay
%     Field
%         :   struct that holds:
%             x
%                 :   The coordinate grid, from -2d to 2d
%             Ex, Ey, Ez, Hx, Hy, Hz
%                 :   All the fields

% default show_plot is false
if nargin < 6
   show_plot = false; 
end

TE = 'TE'; TM = 'TM';

c       = 299792458; % m/s
mu0     = 4*pi*1e-7;
k0      = 2*pi/lambda0;
omega0  = 2*pi*c/lambda0;

% radius of circle u^2 + v^2
V = (pi*d/lambda0)*(nco^2 - ncl^2)^(1/2);

% find # of modes based on radius
n_modes = floor( ( V+(pi/2) ) / (pi/2) );

% Plotting transcendentals
if show_plot
    % coordinates for a circle
    rad = 0:0.01:pi/2;
    x_V = V*cos(rad);
    y_V = V*sin(rad);

    % v = utan(u) and v = -utan(u)
    u1 = 0:0.01:pi/2;
    u2 = pi/2:0.01:pi;
    u3 = pi:0.01:3*pi/2;
    u4 = 3*pi/2:0.01:2*pi;

    v1_te = u1.*tan(u1);
    v2_te = -u2.*cot(u2);
    v3_te = u3.*tan(u3);
    v4_te = -u4.*cot(u4);

    % plotting
    figure;
    h1 = plot(x_V,y_V);
    xlabel('u'); hold on;
    h2 = plot(u1,v1_te,u2,v2_te,u3,v3_te,u4,v4_te);
    ylim([ 0, 2*pi ]);
    xlabel('u'); ylabel('u');
    legend([h2], {'TE even', 'TE odd', 'TE even', 'TE odd'});
    makeFigureNice(); 
end

% Plotting index
if show_plot
    
end

% pre-solve loop init
x_lims  = [0,pi/2];
neff    = [];
beta    = [];
alpha   = [];
kx      = [];

% Initialize field coordinates and variables
x   = linspace(-2*d,2*d,100);
Ex  = zeros(n_modes,length(x));  % field, dimensions of mode x field
Ey  = zeros(n_modes,length(x));
Ez  = zeros(n_modes,length(x));
Hx  = zeros(n_modes,length(x));  % field, dimensions of mode x field
Hy  = zeros(n_modes,length(x));
Hz  = zeros(n_modes,length(x));

%-------------------------------------------------------------------------|
% Solve for the modes
%-------------------------------------------------------------------------|

for ii = 1:n_modes

    if strcmp(pol, TE) % solve for the TE modes

        % find the intersects of the mode/wave equations
        % the intersect is "u"
        if mod(ii,2) == 1 % "even" mode takes on odd mode #
            u = fzero( @(x) symslab_TE_even(x,V), x_lims ); % find the intersects of the mode/wave equations
        else
            u = fzero( @(x) symslab_TE_odd(x,V), x_lims );
        end

    elseif strcmp(pol, TM) % solve for the TM modes

        % find the intersects of the mode/wave equations
        % the intersect is "u"
        if mod(ii,2) == 1 % "even" mode takes on odd mode #
            u = fzero( @(x) symslab_TM_even(x,V,ncl,nco), x_lims ); % find the intersects of the mode/wave equations
        else
            u = fzero( @(x) symslab_TM_odd(x,V,ncl,nco), x_lims );
        end

    end

    % convert "u" to neff and beta
    kx(end+1)   = 2*u/d;  % u = 1/2 kx*d
    beta_ii     = (-kx(ii)^2 + (nco*k0)^2)^(1/2);   % beta
    beta(end+1) = beta_ii;
    neff_ii     = beta_ii/k0;       % neff
    neff(end+1) = neff_ii;
    
    % calculate alpha, A, B, C, D
    if mod(ii,2) == 1 % "even" mode takes on odd mode #
       alpha(end+1) = kx(ii)*tan(kx(ii)*d/2);
       C            = cos(kx(ii)*d/2)/exp(-alpha(end)*d/2);    % asumming A = 1, B = 0
       D            = C;
    else
       alpha(end+1) = -kx(ii)*cot(kx(ii)*d/2); 
       C            = sin(kx(ii)*d/2)/exp(-alpha(end)*d/2);
       D            = -C;
    end
    
    % Solving for field distributions
    if strcmp(pol, TE) % solve for the TE modes
        
        % Make E field
        Ey(ii, x < -d/2)                = C*exp( alpha(ii) * x( x < -d/2) );
        Ey(ii, x >= -d/2 & x <= d/2 )   = cos(kx(ii) * x(x >= -d/2 & x <= d/2) );
        Ey(ii, x > d/2 )                = D*exp( -alpha(ii) * x(x > d/2) );
        
        % Make H field
        Hx(ii, :) = -(beta_ii/(omega0*mu0))*Ey(ii, :);
        
        Hz(ii, x < -d/2)                = (1i/(omega0*mu0)) * alpha(ii) * C * exp( alpha(ii) * x( x < -d/2));
        Hz(ii, x >= -d/2 & x <= d/2)    = -(1i/(omega0*mu0))* kx(ii) * sin( kx(ii) * x(x >= -d/2 & x <= d/2));
        Hz(ii, x > d/2)                 = -(1i/(omega0*mu0)) * alpha(ii) * C * exp( -alpha(ii) * x(x > d/2) );
        
    end
    
    % Plot E field, if plotting is on
%     if show_plot
%         figure;
%         plot( x./d, Ey );
%         xlabel('x (1/d)'); ylabel('E_y(x)');
%         title(sprintf('E_y(x), d = %f \\lambda',d/lambda0));
%         makeFigureNice();
%     end

    x_lims = [ ii*pi/2, (ii+1)*pi/2 ];   % move to next interval

end

Field = struct;
Field.x = x;
Field.Ex = Ex;
Field.Ey = Ey;
Field.Ez = Ez;
Field.Hx = Hx;
Field.Hy = Hy;
Field.Hz = Hz;

end

%-------------------------------------------------------------------------|
% aux functions
%-------------------------------------------------------------------------|

% function to find TE, even modes
function y = symslab_TE_even(u, V)
    y = V^2 - u.^2 - (u.*tan(u)).^2;
end

% function to find TE, odd modes
function y = symslab_TE_odd(u, V)
    y = V^2 - u.^2 - (u.*cot(u)).^2;
end

% function to find TM, even modes
function y = symslab_TM_even(u, V, n1, n2)
    y = V^2 - u.^2 - ( ((n1^2)/(n2^2)) * u.*tan(u)).^2;
end

% function to find zeros of -u*cot(u) = v, and u^2 + v^2 = V^2
function y = symslab_TM_odd(u, V, n1, n2)
    y = V^2 - u.^2 - ( ((n1^2)/(n2^2)) * u.*cot(u)).^2;
end