function [ Xnext ] = transitionFunc( X, dt, params )
    if dt < 0
       Xnext = X;
       return
    end
   
    hexrange = params.footspeed / params.grid.R * dt;
    [si,sj] = nearestHex(X(1,:),X(2,:),params.grid);
    Xnext = zeros(size(X));
    for si_idx=1:size(X,2)
        %[hi,hj] = hexInRange(si(si_idx),sj(si_idx), params.obsmap_hex, hexrange );
        if (si(si_idx) <= 0 || sj(si_idx) <= 0 || si(si_idx) > size(params.hexcache,2) || sj(si_idx) > size(params.hexcache,1) ...
                || isnan(si(si_idx)) || isnan(si(si_idx)) )
            Xnext(:,si_idx) = [nan;nan];
            continue
        end
        
        h = params.hexcache{si(si_idx),sj(si_idx)};

        if ( isempty(h) ) 
            Xnext(:,si_idx) = [nan;nan]; % should kill particle in next observation step
            continue
        end
            
        hi = h(:,1);
        hj = h(:,2);
        
        sel = randi(length(hi));
        Xnext(1,si_idx) = params.grid.xs(hj(sel),hi(sel));
        Xnext(2,si_idx) = params.grid.ys(hj(sel),hi(sel));
    end
    rrand = 0.99*rand(1,size(Xnext,2)) * params.grid.R * sqrt(3)/2;
    trand = rand(1,size(Xnext,2)) * 2*pi;
    Xnext(1,:) = Xnext(1,:) + rrand.*cos(trand);
    Xnext(2,:) = Xnext(2,:) + rrand.*sin(trand);
end