function [ P ] = observeFunc( Y,X,params )
    P = ones(size(X,2),1);
    for i=1:length(params.gaussparams)
        gp = params.gaussparams{i};
        dist = X-repmat(params.rssipos{i}',[1,size(X,2)]);
        dist = sqrt(dist(1,:).^2 + dist(2,:).^2);
        dist = dist / 10;
        muint = interp1(gp(:,1),gp(:,2),dist,'linear','extrap');
        sigint = interp1(gp(:,1),gp(:,2),dist,'nearest','extrap');
        P = P .* normpdf(Y(i),muint,sigint)';
    end
    P(isnan(P)) = 0;
end

