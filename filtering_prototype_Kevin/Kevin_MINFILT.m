function [ AGG_filt2 ] = Kevin_MINFILT( AGG )
%  Implements the MIN filter and successive spike-squashing
N = length(AGG);

MINwin = 6;
AGG_filt = AGG;
for posA=MINwin/2+1:N-MINwin/2
    WL = posA-MINwin/2;
    WR = posA+MINwin/2;
    AGG_filt(posA) = min(AGG(WL:WR));
end



SPK_sizes = {2,4,6,8,10,12,14,16,18,20,22,24,26,28,30};
AGG_filt2 = AGG_filt;
for i=1:length(SPK_sizes)
    AGG_filt_temp = AGG_filt2;
    SPKwin = SPK_sizes{i};
    for posA = SPKwin/2+1:N-SPKwin/2
        WL = posA-SPKwin/2;
        WR = posA+SPKwin/2;
        RNG = AGG_filt2(WL:WR);

        if ( max(RNG) > AGG_filt2(WL) && max(RNG) > AGG_filt2(WR) )
           % flatten the spike
           [zxcv,idx] = max(AGG_filt2(WL:WR));
           AGG_filt_temp(WL:idx+WL) = AGG_filt2(WL);        
           AGG_filt_temp(idx+WL:WR) = AGG_filt2(WR);
        end
        if (min(RNG) < AGG_filt2(WL) && min(RNG) < AGG_filt2(WR))
           % flatten the spike
           [zxcv,idx] = min(AGG_filt2(WL:WR));
           AGG_filt_temp(WL:idx+WL) = AGG_filt2(WL);        
           AGG_filt_temp(idx+WL:WR) = AGG_filt2(WR);
        end
    end
    AGG_filt2 = AGG_filt_temp;
end

end

