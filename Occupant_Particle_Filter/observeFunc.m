function [ P ] = observeFunc( Y,X,params )
    P = ones(size(X,2),1);
    for i=1:length(params.gaussparams)
        if ~isnan(Y(i))
            gp = params.gaussparams{i};
            dist = X-repmat(params.rssipos{i}',[1,size(X,2)]);
            dist = sqrt(dist(1,:).^2 + dist(2,:).^2);
            dist = dist / 10;
            muint = interp1(gp(:,1),gp(:,2),dist,'linear','extrap');
            %sigint = interp1(gp(:,1),gp(:,3),dist,'nearest','extrap');
            %muint = muint + 20;
            sigint = 40;
            % gaussian prob
            P = P .* normpdf(Y(i),muint,sigint)';
            
            % uniform prob
            %P = P .* (( (Y(i) >= muint-sigint) & (Y(i) <= muint + sigint) )./(2*sigint))';
            
        end
    end
    P(isnan(P)) = 0;
end