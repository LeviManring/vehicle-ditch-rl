function bounds = bounds_fun(rows)
%% bounds_fun
% Levi Manring, Duke University
% 2021
%
% This function calculates the bounds on a region, used in
% bestresultsplotting.m
%
% Inputs:
%   rows: nx1 array indicating the rows of an array that meet some
%   condition
%
% Outputs:
%   bounds: mx2 array indicating the boundaries of the condition. Each row
%           contians 2 values, the start of the boundary and the end of the
%           boundary

bounds = rows(1,1);
bounds_k = 1;
column = 2;
for k = 1:(length(rows) - 1)
    if rows(k+1,1) > rows(k,1) + 1
        bounds(bounds_k,column) = rows(k,1);
        column = column - 1;
        bounds(bounds_k + 1,column) = rows(k+1,1);
        bounds_k = bounds_k + 1;
        column = column + 1;
    end
end
bounds(bounds_k,column) = rows(end,1);