% testing AIM grating cell

clear; close all;

dxy     = 10;
lambda  = 1550;
k0      = 2*pi/lambda;

neff_slab   = 2.842431701;
n_clad      = 1.45;
theta       = 15; % deg
period      = 2*pi/( k0 * ( neff_slab - n_clad*sin( theta*pi/180 ) ) );
period      = round(period/dxy)  * dxy;

fill_top = 0.9;
fill_bot = 0.9;

offset_ratios = linspace(0,0.99, 10);

% saving variables
directionality_vs_offset = zeros(size(offset_ratios));

% pick which grating design to try
which_grating = 'default si custom sin gratings full etch';

% set grating solver settings
num_modes   = 15;
BC          = 0;    % 0 = PEC
pml_options = [1, 100, 20, 2]; 
OPTS        = struct();
guessk      = k0 * neff_slab;

tic;
for i_offset = 1:length(offset_ratios)
    
    fprintf('loop %i of %i\n', i_offset, length(offset_ratios));

    switch which_grating
        case 'default si custom sin gratings full etch'
            % custom version that ovverides the sin thickness and spacings
            % but assumes fully etched sin layers (dual layers)
            % fields will be "si_to_bot_sin", "bot_sin_thick",
            % "sin_to_sin_thick", "top_sin_thick"

            si_to_bot_sin       = 100;
            bot_sin_thick       = 60;
            sin_to_sin_thick    = 0;
            top_sin_thick       = 160;

            OPTS = struct( 'geometry','default si custom sin gratings full etch', ...
                           'si_to_bot_sin', si_to_bot_sin, ...
                           'bot_sin_thick', bot_sin_thick, ...
                           'sin_to_sin_thick', sin_to_sin_thick, ...
                           'top_sin_thick', top_sin_thick );

            GC = f_makeGratingCell_AIM( dxy, lambda, ...
                                    period, fill_top, fill_bot, offset_ratios(i_offset), OPTS  );
    end

    % run simulation
    GC = GC.runSimulation( num_modes, BC, pml_options, guessk, OPTS );
    
    % update k
    guessk = GC.k;
    
    % save directionality
    directionality_vs_offset(i_offset) = GC.directivity;
    
    toc;
    
end

% plot e field
% GC.plot_E_field_gui();

% directionality vs. offset
figure;
plot( offset_ratios, directionality_vs_offset, '-o' );
xlabel('offset ratio'); ylabel('directionality up/down');
makeFigureNice();

















