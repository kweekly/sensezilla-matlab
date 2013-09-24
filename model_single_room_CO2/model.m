
[tout,yout] = ode45(@transfun,[0 100],10);
plot(tout,yout);