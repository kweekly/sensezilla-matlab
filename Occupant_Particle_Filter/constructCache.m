function [ rangecache ] = constructCache( obsmap, range )
%  [ rangecache ] = constructCache( obsmap, range )
%   Generate a cache of hexes in range
%   i.e. rangecache{i,j} = [[1 2];[3 4]...] giving list of reachable i,j
%   pairs

rangecache = cell(size(obsmap));
for i=1:size(obsmap,2)
    fprintf(['i=' num2str(i) '/' num2str(size(obsmap,2)) '\n']);
    for j=1:size(obsmap,1)
        [hi,hj] = hexInRange( i, j, obsmap, range );
        rangecache{i,j} = [hi,hj];
    end
end
end

