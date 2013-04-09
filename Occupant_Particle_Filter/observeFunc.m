function [ P ] = observeFunc( Y,X,params )
    P = (normpdf(Y(1),X(1,:),params.wifiSigma) .* normpdf(Y(2),X(2,:),params.wifiSigma))';
    P(isnan(P)) = 0;
    P(isnan(Y)) = 1; % for missing data
end

