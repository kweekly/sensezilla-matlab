%%


imgbounds = [[-160.11,468.5];[-308.07,159.54]];
[obsmap, entrymap, stopmap] = readFloorplan('floorplan.png');
grid = hexGrid(imgbounds(1,:),imgbounds(2,:),2,4);
obsmap_hex = rect2hexgrid(obsmap, grid);
entrymap_hex = rect2hexgrid(entrymap, grid);
stopmap_hex = rect2hexgrid(stopmap, grid);

%%
%{
params.grid = grid;
params.entrymap_hex = entrymap_hex;
params.stopmap_hex = stopmap_hex;
params.obsmap_hex = obsmap_hex;
params.wifiSigma = 50;
params.footspeed = 1.5;
params.dt = 60; 
params.hexspeed = params.footspeed / params.grid.R * params.dt;

params.hexcache = constructCache(obsmap_hex, params.hexspeed);
save('params.mat','params');
%}
load('params.mat');

%%
dat = importdata('data/test.csv');
data = dat.data(dat.data(:,2)==1,:);


cut = 10;
data = data(1:cut,:);

data(:,4) = -data(:,4);
Tdata = data(:,1);
Ydata = [data(:,3) data(:,4)]';

% add in extraneous data
%{
Tdt = 10;
t_temp = t;
Y_temp = Y;
t = [min(t):Tdt:max(t)]';
Y = nan*ones(2,size(t,1));
ttidx = 1;
for tidx=1:length(t)-1
    if t_temp(ttidx) >= t(tidx) && t_temp(ttidx) < t(tidx + 1)
        Y(:,tidx) = Y_temp(:,ttidx);
        ttidx = ttidx + 1;
    end
end
%}

%%

f = figure;
set(f,'Renderer','zbuffer');
hold on;
[obsj,obsi] = find(obsmap_hex > 0.2);
hplot = hexPlot(obsi,obsj,'k',grid);
%surf(grid.xs,grid.ys,entrymap_hex-3,'EdgeColor','None');
%surf(grid.xs,grid.ys,stopmap_hex-1,'EdgeColor','None');
view([0, 90]);
visited = unique(data(:,[3 4]),'rows');
[visi,visj] = nearestHex(visited(:,1),visited(:,2),grid);
hexPlot(visi,visj,'r',grid); axis equal;

%% Uniform Stupid Sampling

% Create uniform samples
nSamples = 10000;
Xest = rand([2 nSamples]);
Xest(1,:) = Xest(1,:) * (imgbounds(1,2) - imgbounds(1,1)) + imgbounds(1,1);
Xest(2,:) = Xest(2,:) * (imgbounds(2,2) - imgbounds(2,1)) + imgbounds(2,1);

West = ones(nSamples,1) * 1/nSamples;


T = unique([min(Tdata):params.dt:max(Tdata) max(Tdata)]');
Y = interp1(Tdata,Ydata',T,'linear')';


% Transition function
%transfunc = @(X,dt) X + dt*wSpeed*2*(rand(size(X))-0.5);
% Observation function
%observefunc = @(Y,X) normpdf(Y(1),X(1,:),wifiSigma).*normpdf(Y(2),X(2,:),wifiSigma);
% do the filtering

[Xest,West,Xtraj] = SIRFilter(Xest,West,Y,T, @transitionFunc, @observeFunc, params);


close all;

figure;
subplot(2,1,1); hold on; plot(T,Y(1,:),'+-b'); plot(T,Xtraj(1,:),'o-g');
subplot(2,1,2); hold on; plot(T,Y(2,:),'+-b'); plot(T,Xtraj(2,:),'o-g');


f = figure;
set(f,'Renderer','zbuffer');
hold on;

[obsj,obsi] = find(obsmap_hex > 0.2);
%{
for ii=1:length(obsj)
   plot(grid.xs(obsj(ii),obsi(ii)),grid.ys(obsj(ii),obsi(ii)),'r.','MarkerSize',4); 
end
%}
hplot = hexPlot(obsi,obsj,'k',grid);
view([0, 90]);
axis equal;


plot(Xest(1,:),Xest(2,:),'g.','MarkerSize',1);

plot(Y(1,:),Y(2,:),'ro-');

%visited = unique(data(:,[3 4]),'rows');
%[visi,visj] = nearestHex(visited(:,1),visited(:,2),grid);
%hexPlot(visi,visj,'r',grid); axis equal;

%[esti,estj] = nearestHex(Xtraj(1,:)',Xtraj(2,:)',grid);
%hexPlot(esti,estj,'b',grid); axis equal;
plot(Xtraj(1,:),Xtraj(2,:),'b-x');


