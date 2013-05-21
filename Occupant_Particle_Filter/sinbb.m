%%


imgbounds = [[-160.11,468.5];[-308.07,159.54]];
[obsmap, entrymap, stopmap] = readFloorplan('floorplan.png');
grid = hexGrid(imgbounds(1,:),imgbounds(2,:),2,4);
obsmap_hex = rect2hexgrid(obsmap, grid);
entrymap_hex = rect2hexgrid(entrymap, grid);
stopmap_hex = rect2hexgrid(stopmap, grid);

%% Load params and construct map

params.grid = grid;
params.entrymap_hex = entrymap_hex;
params.stopmap_hex = stopmap_hex;
params.obsmap_hex = obsmap_hex;
params.wifiSigma = 50;
params.footspeed = 1.5;
params.dt = 5; 
params.hexspeed = params.footspeed / params.grid.R * params.dt;

params.hexcache = constructCache(obsmap_hex, params.hexspeed);
save('params.mat','params');

load('params.mat');

%% Load sensor data
rssi_origin = [20 -132];
cd data/RSSIdata
rssi_system;

gaussparams = {};
rssidata = {};
rssipos = {};

keys = sense_pos.keys;
for ki=1:length(keys)
    senseid = keys{ki};
    rssipos{ki} = sense_pos(senseid);
    gaussparams{ki} = csvread(['gaussparams_s' senseid '.csv']);
    dat = importdata(['RSSI_t0005_s' senseid '.csv']);
    if ( isstruct(dat) )
        rssidata{ki} = dat.data;
    else
       rssidata{ki} = zeros(0,2); 
    end
end
cd ../..

params.rssipos = rssipos;
params.gaussparams = gaussparams;

T = [];
for ki=1:length(keys)
   d = rssidata{ki};
   T = union(T, d(:,1));
end
T = T / 1000;

Y = nan*ones(length(keys),length(T));
for ki=1:length(keys)
    d = rssidata{ki};
    for ti=1:size(d,1)
        Y(ki,T==d(ti,1)) = d(ti,2);
    end
end

%% Load Ground truth data
gtdata = importdata('data/ground_truth.txt');
gtdata_t = gtdata(:,1) + T(1);
gtdata_x = [rssi_origin(1)+gtdata(:,3)*10, rssi_origin(2) + gtdata(:,2)*10];

%% Uniform Stupid Sampling

% Create uniform samples
nSamples = 1000;
Xest = rand([2 nSamples]);
Xest(1,:) = Xest(1,:) * (imgbounds(1,2) - imgbounds(1,1)) + imgbounds(1,1);
Xest(2,:) = Xest(2,:) * (imgbounds(2,2) - imgbounds(2,1)) + imgbounds(2,1);

West = ones(nSamples,1) * 1/nSamples;

% Transition function
%transfunc = @(X,dt) X + dt*wSpeed*2*(rand(size(X))-0.5);
% Observation function
%observefunc = @(Y,X) normpdf(Y(1),X(1,:),wifiSigma).*normpdf(Y(2),X(2,:),wifiSigma);
% do the filtering

[Xest,West,Xtraj] = SIRFilter(Xest,West,Y,T, @transitionFunc, @observeFunc, params);


close all;

figure;
subplot(2,1,1); hold on; plot(gtdata_t,gtdata_x(:,1),'+-b'); plot(T,Xtraj(1,:),'o-g');
subplot(2,1,2); hold on; plot(gtdata_t,gtdata_x(:,2),'+-b'); plot(T,Xtraj(2,:),'o-g');


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
plot(gtdata_x(:,1),gtdata_x(:,2),'LineWidth',2);

plot(Xest(1,:),Xest(2,:),'g.','MarkerSize',1);

%visited = unique(data(:,[3 4]),'rows');
%[visi,visj] = nearestHex(visited(:,1),visited(:,2),grid);
%hexPlot(visi,visj,'r',grid); axis equal;

%[esti,estj] = nearestHex(Xtraj(1,:)',Xtraj(2,:)',grid);
%hexPlot(esti,estj,'b',grid); axis equal;
plot(Xtraj(1,:),Xtraj(2,:),'b-x');


