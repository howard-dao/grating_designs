function [ poly_min_len, poly_min_gap, ...
           body_min_len, body_min_gap, ...
           par_etch_min_len, par_etch_min_gap, ...
           sin_min_len, sin_min_gap ] = f_clo_min_feat_sizes()
% Returns min. feature size for clo process
%
% all units in nm

poly_min_len = 100;
poly_min_gap = 100;

% full etch
body_min_len = 80;
body_min_gap = 90;

% partial etch
par_etch_min_len = 150; % defined as min KG width
par_etch_min_gap = 200; % defined as min KG space

% nitride (potentially not updated for v1.0.2.0)
sin_min_len = 120;
sin_min_gap = 190;


end

