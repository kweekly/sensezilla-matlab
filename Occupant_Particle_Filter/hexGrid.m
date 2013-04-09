function [ grid ] = hexGrid( X, Y, R, COV )
%HEXGRID create hexagonal grid structure
%  grid = hexGrid(X, Y, R) returns a structure
%  representing a hexagonally-tiled grid, where the hexagons have radius R.
%  X and Y are two-element vectors giving the bounds that should be tiled.
%  The following fields are in the returned structure:
%           grid.size : A two element vector giving number of entries of
%           the hex tiles in the x and y direction
%           grid.xs : x-coordinate of centerpoint of hex tiles
%           grid.ys : y-coordinate of centerpoint of hex tiles
%           grid.R  : the radius of the hexagons ( input argument )
%           grid.rectbounds : The original input arguments
%           [X(1),X(2);Y(1),Y(2)]
%           grid.covbounds : The covered area
%
% grid = hexGrid(X,Y,R,COV)  indicates that "(COV)" extra space should be
% given to overlap region
%
% See: http://blog.ruslans.com/2011/02/hexagonal-grid-math.html

grid.rectbounds = [[X(1) X(2)];[Y(1) Y(2)]];
grid.R = R;

if nargin == 4 && COV > 0
    X(1) = X(1) - COV*R/2;
    X(2) = X(2) + COV*R/2;
    Y(1) = Y(1) - COV*R/2;
    Y(2) = Y(2) + COV*R/2;
end

grid.covbounds = [[X(1) X(2)];[Y(1) Y(2)]];

H = 2*R*sin(60*pi/180);
S = 3/2*R;
W = 2*R;

nTiles_x = floor(((X(2) - X(1)) - W)/S) + 1;
nTiles_y = floor(((Y(2) - Y(1)) - H/2)/H) + 1;
if (nTiles_x <= 0 || nTiles_y <= 0)
   fprintf(2, ['Error: Hex tiles are too big to fit in given x/y bounds\n']); 
   return;
end

grid.size = [nTiles_x;nTiles_y];
[grid.xs,grid.ys] = meshgrid(0:nTiles_x-1,0:nTiles_y-1);
grid.ys = grid.ys * H - mod(grid.xs,2) * H/2 - H/2 + Y(1);
grid.xs = grid.xs * S + R + X(1);
end