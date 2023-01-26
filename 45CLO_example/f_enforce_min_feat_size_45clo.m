function is_ok = f_enforce_min_feat_size_45clo( period, top_fill, bot_fill, offset, etch_depth, ...
                                                custom_pc_minlen )
% Returns true if datapoint satisfies min. feature size demands, otherwise
% returns false if it violates
%
% User can add arguments that determine the min. feature sizes
%
% Inputs:
%   period
%       type: double, scalar
%       desc: period of unit cell
%   top_fill
%       type: double, scalar
%       desc: ratio of top layer length/period
%   bot_fill
%       type: double, scalar
%       desc: ratio of bot layer length/period
%   offset
%       type: double, scalar
%       desc: offset of bottom layer from left edge of cell
%   etch_depth
%       type: string
%       desc: either 'shallow' or 'full'
%   custom_pc_minlen
%       type: scalar
%       desc: custom pc min. length

is_ok = true;

% get min feat sizes for clo
[ poly_min_len, poly_min_gap, ...
   body_min_len, body_min_gap, ...
   par_etch_min_len, par_etch_min_gap, ...
   sin_min_len, sin_min_gap ] = f_clo_min_feat_sizes();

% override pc min length
if nargin == 6
    poly_min_len = custom_pc_minlen;
end

switch etch_depth
    case 'shallow'
        bot_min_len = par_etch_min_gap;
        bot_min_gap = par_etch_min_len;
    case 'full'
        bot_min_len = body_min_len;
        bot_min_gap = body_min_gap;
end
 

% check for min. waveguide length is satisfied
if period * top_fill < poly_min_len || period * bot_fill < bot_min_len
    is_ok = false;
end

% check for min. gap is satisfied
if period * (1-top_fill) < poly_min_gap || period * (1-bot_fill) < bot_min_gap
    is_ok = false;
end
            
end
