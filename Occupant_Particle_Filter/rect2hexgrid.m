function [ H ] = rect2hexgrid( X,Y,A, hexgrid )
%RECT2HEXGRID 
%   [ H ] = rect2hexgrid( X,Y,A, hexgrid )
%Interpolate a matrix in a certain rectangular coordinate
%system to the hex grid specified. If X and Y are vectors, they will be
%'meshgrid'ed.
%
% If X and Y are not specified, then bounds are read from hexgrid.

if (nargin == 2)
   A = X;
   hexgrid = Y;
   X = hexgrid.rectbounds(1,1):((hexgrid.rectbounds(1,2)-hexgrid.rectbounds(1,1))/size(A,2)):hexgrid.rectbounds(1,2);
   X = X(1:end-1);
   Y = hexgrid.rectbounds(2,1):((hexgrid.rectbounds(2,2)-hexgrid.rectbounds(2,1))/size(A,1)):hexgrid.rectbounds(2,2);
   Y = Y(1:end-1);
end

if (isvector(X) && isvector(Y))
   [X,Y] = meshgrid(X,Y); 
end

if ( all(size(X) ~= size(A)) || all(size(Y) ~= size(A)) )
   error('Size of X and Y and A must match');
end

H = interp2(X,Y,A,hexgrid.xs,hexgrid.ys);
end

