function [ Xnext ] = transitionFunc( X, dt, params )
    hexrange = params.footspeed / params.grid.R * dt;
    [si,sj] = nearestHex(X(1,:),X(2,:),params.grid);
    Xnext = zeros(size(X));
    for si_idx=1:size(X,2)
        [hi,hj] = hexInRange(si(si_idx),sj(si_idx), params.obsmap_hex, hexrange );
        if ( isempty(hi) ) 
            Xnext(:,si_idx) = [nan;nan]; % should kill particle in next observation step
            continue
        end
            
        sel = randi(length(hi));
        Xnext(1,si_idx) = params.grid.xs(hj(sel),hi(sel));
        Xnext(2,si_idx) = params.grid.ys(hj(sel),hi(sel));
    end
end