function [dataf] = do_calibrate_sametime(time, data, ctime, cdata)
    dataf = zeros( [length(ctime) 1] );
    j = 1;
    for i=1:length(ctime)-1
        j2 = find(time(j+1:end) >= ctime(i+1),1,'first') + j;
        dataf(i) = mean(data(j:j2));
        j = j2;
    end
    dataf(end) = dataf(end-1);
end