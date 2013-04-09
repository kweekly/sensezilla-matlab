function [ i,j ] = nearestHex( x,y, grid )
%[i,j] = NEARESTHEX(x,y,grid)
%  Find the hex indices for x and y points

R = grid.R;
H = 2*R*sin(60*pi/180);
S = 3/2*R;
W = 2*R;

i = floor((x - grid.covbounds(1,1))/S) + 1;
j = ceil(((y - grid.covbounds(2,1)) + mod(i,2)*H/2 - H/2)/H) + 1;

end

