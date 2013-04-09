function [datf] = do_calibrate_filter(dat, nsamp)
    datf = conv(dat, ones([nsamp 1])./nsamp, 'same');
end