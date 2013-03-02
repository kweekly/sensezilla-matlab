function [ Xnext ] = transitionFunc( X, dt )
    wSpeed = 10;
    R = wSpeed * rand(1,size(X,2)) * dt;
    T = 2*pi*rand(1,size(X,2));
    Xnext = X + [R.*cos(T);R.*sin(T)];
end

