function h = plot_bz_slice(bzSlice, plot_args)
% h = plot_bz_slice(bzSlice, plot_args)
%   This is a function that plots the planar Brilluoin zone 
%   extracted from the 'get_bz_slice()' function.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   bzSlice:        MATLAB data structure containing BZ slice data from 'get_bz_slice()'
%   -   plot_args:      plot arguments for the Brilluoin zone slice
%
%   OUT: (none)

%% Default parameters
if nargin < 2; plot_args = []; end
if isempty(plot_args); plot_args = []; end
%% - 1 - Initialising the transformation parameters
hold on;
ax = gca;
XLim = ax.XLim; YLim = ax.YLim; 
for i = 1:length(bzSlice.X)
    h(i) = plot(bzSlice.X{i}, bzSlice.Y{i}, 'k-');
end
axis([XLim, YLim]);
end