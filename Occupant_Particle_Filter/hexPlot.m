function [ p ] = hexPlot( i, j, color, grid )
%p = HEXPLOT(i,j, color, grid)
%   fill in hexes on the grid using "patch" objects
%   i and j are indices of the hexagons on the "grid" argument
%   returns handle of patch created

R = grid.R;
H = 2*R*sin(60*pi/180);
S = 3/2*R;
W = 2*R;


centx = (i-1) * S + R + grid.covbounds(1,1);
centy = (j-1) * H - mod(i,2) * H/2 + grid.covbounds(2,1);
deg = [0:pi/3:2*pi];
xverts = repmat(centx,1,7) + repmat(R*cos(deg),length(centx),1);
yverts = repmat(centy,1,7) + repmat(R*sin(deg),length(centx),1);
p = patch(xverts',yverts',color);
end

