%%


imgbounds = [[-160.11,468.5];[-308.07,159.54]];
[obsmap, entrymap, stopmap] = readFloorplan('floorplan.png');
grid = hexGrid(imgbounds(1,:),imgbounds(2,:),2,4);
obsmap_hex = rect2hexgrid(obsmap, grid);
entrymap_hex = rect2hexgrid(entrymap, grid);
stopmap_hex = rect2hexgrid(stopmap, grid);

%%
dat = importdata('data/test.csv');
data = dat.data(dat.data(:,2)==1,:);


cut = 10;
data = data(1:cut,:);

data(:,4) = -data(:,4);
t = data(:,1);
Y = [data(:,3) data(:,4)]';

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

params.grid = grid;
params.entrymap_hex = entrymap_hex;
params.stopmap_hex = stopmap_hex;
params.obsmap_hex = obsmap_hex;
params.wifiSigma = 250;
params.footspeed = 0.5;
params.dt = 1;

% Transition function
%transfunc = @(X,dt) X + dt*wSpeed*2*(rand(size(X))-0.5);
% Observation function
%observefunc = @(Y,X) normpdf(Y(1),X(1,:),wifiSigma).*normpdf(Y(2),X(2,:),wifiSigma);
% do the filtering
[Xest,West,Xtraj] = SIRFilter(Xest,West,Y,t, @transitionFunc, @observeFunc, params);
figure;
subplot(2,1,1); hold on; plot(t,Y(1,:),'+-b'); plot(t,Xtraj(1,:),'o-g');
subplot(2,1,2); hold on; plot(t,Y(2,:),'+-b'); plot(t,Xtraj(2,:),'o-g');

%%

f = figure;
set(f,'Renderer','zbuffer');
hold on;
[obsj,obsi] = find(obsmap_hex > 0.2);
hplot = hexPlot(obsi,obsj,'k',grid);
view([0, 90]);
plot(Y(1,:),Y(2,:),'ro-');

%visited = unique(data(:,[3 4]),'rows');
%[visi,visj] = nearestHex(visited(:,1),visited(:,2),grid);
%hexPlot(visi,visj,'r',grid); axis equal;

[esti,estj] = nearestHex(Xtraj(1,:)',Xtraj(2,:)',grid);
hexPlot(esti,estj,'b',grid); axis equal;
plot(Xtraj(1,:),Xtraj(2,:),'b-');


