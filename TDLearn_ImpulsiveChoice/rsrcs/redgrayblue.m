function c = redgrayblue(m)
%REDGRAYBLUE    Shades of red and blue color map
%   REDGRAYBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to gray, and then through shades of red to bright red.
%   REDGRAYBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redgrayblue)
%
%   See also REDBLUE, HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   EHS:20160107

if nargin < 1, m = size(get(gcf,'colormap'),1); end

grayVal = 0.5;

if (mod(m,2) == 0)
    % From [0 0 1] to [0.5 0.5 0.5], then [0.5 0.5 0.5] to [1 0 0];
    m1 = m*0.5;
    r = ((0:m1-1)'/max(m1-1,1)).*grayVal;
    g = r;
    r = [r; r+ones(m1,1).*(1-grayVal)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [0.5 0.5 0.5], then [0.5 0.5 0.5] to [1 0 0];
    m1 = floor(m*0.5);
    r = ((0:m1-1)'/max(m1,1)).*grayVal;
    g = r;
    r = [r(1:end-1); linspace(max(r),1,m1+2)'];
    g = [g; grayVal; flipud(g)];
    b = flipud(r);
end

c = [r g b]; 
