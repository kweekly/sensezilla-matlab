function [ P ] = observeFunc( Y,X )
    wifiSigma = 1000;
    P = (normpdf(Y(1),X(1,:),wifiSigma) .* normpdf(Y(2),X(2,:),wifiSigma))';
end

