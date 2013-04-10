function [ endi,endj ] = hexInRange( starti, startj, obsmap, range )
%[endi,endj] = HEXINRANGE( starti, startj, obsmap, range )
%   Return set of hexes within range of hex at starti,startj while avoid
%   obstacles in obsmap

neighbor0 = [[-1,-1];[0,-1];[1,-1];[1,0];[0,1];[-1,0]];
neighbor1 = [[0,1];[0,0];[0,1];[0,1];[0,0];[0,1]];

if ( starti > size(obsmap,2) || starti <= 0 || startj > size(obsmap,1) || startj <= 0 || obsmap(startj,starti) > 0.2) 
    endi = [];
    endj = [];
    return
end
newset = [starti startj];
retval = zeros(0,2);
for i=1:range
    candidates = zeros(0,2) ;
    for j=1:size(newset,1)
        candidates = union(candidates,repmat(newset(j,:),size(neighbor0,1),1) + neighbor0 + mod(newset(j,1),2)*neighbor1,'rows');
    end
    candidates = setdiff(candidates,union(retval,newset,'rows'),'rows');
    
    remval = zeros(0,2);
    
    for k=1:size(candidates,1)
       if (any(candidates(k,:) <= 0) || candidates(k,2) > size(obsmap,1) || candidates(k,1) > size(obsmap,2) || obsmap(candidates(k,2),candidates(k,1)) > 0.2)
           remval = [remval;candidates(k,:)];
       end
    end
    
    retval = union(retval,newset,'rows');
    newset = setdiff(candidates,remval,'rows');
end

retval = union(retval,newset,'rows');
endi = retval(:,1);
endj = retval(:,2);
end