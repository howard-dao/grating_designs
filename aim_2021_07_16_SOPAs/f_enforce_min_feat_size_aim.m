function is_ok = f_enforce_min_feat_size_aim(   period, top_fill, ...
                                                bot_fill, offset ...
                                                )
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

is_ok = true;

[   bot_sin_min_len, bot_sin_min_gap, ...
    top_sin_min_len, top_sin_min_gap ] = f_aim_min_feat_sizes();

% check if minimum waveguide length is satisfied
if period * top_fill < top_sin_min_len || period * bot_fill < bot_sin_min_len
    is_ok = false;
end

% check if minimum gap is satisfied
if period * (1-top_fill) < top_sin_min_gap || period * (1-bot_fill) < bot_sin_min_gap
    is_ok = false;
end

end
