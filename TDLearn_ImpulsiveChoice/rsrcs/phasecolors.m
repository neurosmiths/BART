function c = phasecolors()
%PHASECOLORS Shades of red and blue color map
%   PHASECOLORS(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with black and range through shades of
%   blue to white, and then through shades of red back to black.
%   PHASECOLORS, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(phasecolors)
%
%   See also REDBLUE, HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   EHS:20160107

if nargin < 1, m = size(get(gcf,'colormap'),1); end

load('/home/user1/code/matlab/analyzeMSIT/dependencies/phaseColormap.mat')
c = phaseCmap;
