function [Z, counts, xbins, ybins] = histheatmap(x,y,xbins,ybins)
    if isscalar(xbins)
       xbins = min(x):(max(x)-min(x))/xbins:max(x); 
    end

    if isscalar(ybins)
       ybins = min(y):(max(y)-min(y))/ybins:max(y); 
    end

    %[X, Y] = meshgrid(xbins,ybins);
    Z = zeros([length(ybins) length(xbins)]);
    counts = zeros([1 length(xbins)]);
    
    for i=1:length(xbins)-1
       yfilt = y(x >= xbins(i) & x < xbins(i+1));
       Z(:,i) = histc(yfilt,ybins);
       counts(i) = sum(Z(:,i));
       if (counts(i) > 0)
         %Z(:,i) = Z(:,i) / counts(i);
       else
         Z(:,i) = 0;
       end
    end
    
    cmap = colormap;
    mmax = max(max(Z));
    for i=1:length(xbins)-1
        mmax = max(1e-25,max(Z(:,i)));        
        for j=1:length(ybins)-1
            x2 = xbins(i+1);
            y2 = ybins(j+1);
            if ( isinf(x2) )
                x2 = xbins(i) + (xbins(i)-xbins(i-1));
            end
            if ( isinf(y2) )
                y2 = ybins(j) + (ybins(j)-ybins(j-1));
            end
            patch([xbins(i) x2 x2 xbins(i)],...
                   [ybins(j) ybins(j) y2 y2],...
                   [Z(j,i) Z(j,i) Z(j,i) Z(j,i)],...
                   cmap(1+floor( Z(j,i)/mmax * (size(cmap,1)-1)),:));
        end
    end
    
    view([0 90]);
    %m = surf(X, Y, Z, 'EdgeColor','none');
    %view([0 90]);
end