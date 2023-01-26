function [] = f_design_45CLO_gc_duallayer(filepath, fill_top_override, fill_bot_override, ...
                                                            fill_vs_top_bot_both, etch_depth, MFD, apodized, ...
                                                            run_fdtd, fdtd_opts )
% synthesizes final design
%
% inputs
%   filepath
%       path to where synth_obj_45CLO.mat is
%   fill_top_override
%       for unit cell quadrant selection
%   fill_bot_override
%       for unit cell quadrant selection
%   fill_vs_top_bot_both
%       for unit cell selection. either 'top', 'bottom', or 'both'
%   etch_depth
%       determines RX etch depth. either 'shallow' or 'full'
%   MFD
%       units nm
%   apodized
%       if true, design apodized, otherwise uniform
%   run_fdtd
%       if true, run fdtd and save coupling results
%   fdtd_opts
%       struct, see f_run_fdtd_2Dgrating_clo

filename = 'synth_obj_45CLO.mat';

% input wg options
input_wg_type           = 'bottom'; % for synthesis object
wg_interface_type       = 'body'; % for fdtd

% load synth_obj
synthdata = load( [ filepath filesep filename ] );
synth_obj = synthdata.synth_obj;

% function to enforce min feat sizes
enforce_min_feat_size_func = @( period, top_fill, bot_fill, offset ) ...
 f_enforce_min_feat_size_45clo( period, top_fill, bot_fill, offset, etch_depth, 100 );

% generate design
if apodized
    synth_obj = synth_obj.generate_final_design_apodized_gaussian( MFD, ...
                                input_wg_type, ...
                                enforce_min_feat_size_func, ...
                                fill_top_override, ...
                                fill_bot_override, ...
                                fill_vs_top_bot_both );
    design_label = [ datestr(now,'yymmdd_HHMM') '_apod_mfd' num2str(MFD*1e-3) ];
else
    synth_obj = synth_obj.generate_final_design_uniform_gaussian( MFD, ...
                                input_wg_type, ...
                                enforce_min_feat_size_func, ...
                                fill_top_override, ...
                                fill_bot_override );
    design_label = [ datestr(now,'yymmdd_HHMM') '_unif_mfd' num2str(MFD*1e-3) ];
end

% plot final design
f_plot_final_design( synth_obj, 'dual' );

% make path to save final design to
save_data_path = [ filepath filesep 'final_design' filesep design_label ];
mkdir( save_data_path );
                                                        
% save the final synthesis object and workspace
save( [ save_data_path filesep 'synthobjfinal' ] );

if run_fdtd
    
    % run fdtd
    tic;
    [ fdtd_results ] = f_run_fdtd_2Dgrating_clo( synth_obj, etch_depth, ...
                                            'dual', wg_interface_type, [ save_data_path filesep '2dfdtd.fsp' ], fdtd_opts );
    toc;

    % calculate coupling results
    [ coupling_results ] = f_calc_coupling_results( synth_obj, fdtd_results );

    % save data
    save( [ save_data_path filesep 'final_results.mat' ] );

    % plot results
    f_plot_final_results(save_data_path, 3);

    save_all_figs( save_data_path );
    
end

end
























