function [ O, E, S ] = readFloorplan( file_name )
%READFLOORPLAN Read floorplan file in png into matrices
%       [ O, E, S ] = readFloorplan( file_name ) : reads from the png file
%       given by file_name into O (obstacles), E ( entry/exits ), S (
%       static areas, i.e. places where someone would stay still )

rgb = imread(file_name);
O = (rgb(:,:,1) == rgb(:,:,2)) & (rgb(:,:,1) == rgb(:,:,3)) & (rgb(:,:,1) ~= 255);
E = (rgb(:,:,1) == 0) & (rgb(:,:,2) > 10) & (rgb(:,:,3) == 0);
S = (rgb(:,:,1) > 10) & (rgb(:,:,2) > 10) & (rgb(:,:,3) == 0);
O = O(end:-1:1,1:end);
E = E(end:-1:1,1:end);
S = S(end:-1:1,1:end);
end

