%%


imgbounds = [[-160.11,468.5];[-308.07,159.54]];
[obsmap, entrymap, stopmap] = readFloorplan('floorplan.png');
grid = hexGrid(imgbounds(1,:),imgbounds(2,:),2,4);
obsmap_hex = rect2hexgrid(obsmap, grid);
entrymap_hex = rect2hexgrid(entrymap, grid);
stopmap_hex = rect2hexgrid(stopmap, grid);

%% Load params and construct map
%{
params.grid = grid;
params.entrymap_hex = entrymap_hex;
params.stopmap_hex = stopmap_hex;
params.obsmap_hex = obsmap_hex;
params.wifiSigma = 50;
params.footspeed = 1.5;
params.dt = 5; 
params.hexspeed = params.footspeed / params.grid.R * params.dt;

%params.hexcache = constructCache(obsmap_hex, params.hexspeed);
save('params.mat','params');
%}
load('params.mat');

%% Load Ground truth data
rssi_origin = [20 -132];
gtdata = importdata('data/ground_truth.txt');
gtdata_t = gtdata(:,1);
gtdata_x = [rssi_origin(1)+gtdata(:,3)*10, rssi_origin(2) + gtdata(:,2)*10];

%% Load sensor data

cd data/RSSIdata
rssi_system;

gaussparams = {};
rssidata = {};
rssipos = {};

keys = sense_pos.keys;
mint = 9e99;
for ki=1:length(keys)
    senseid = keys{ki};
    sp = sense_pos(senseid);
    rssipos{ki} = [rssi_origin(1)+sp(2)*10, rssi_origin(2) + sp(1)*10];
    gp = csvread(['gaussparams_s' senseid '.csv']);
    dat = importdata(['RSSI_t0006_s' senseid '.csv']);
    
    gaussparams{ki} = gp;
    
    if ( isstruct(dat) )
        rssidata{ki} = dat.data;
        if min(dat.data(:,1))<mint
            mint = min(dat.data(:,1));
        end
    else
       rssidata{ki} = zeros(0,2); 
    end
end
cd ../..

for ki=1:length(keys)
   rssidata{ki}(:,1) = rssidata{ki}(:,1)-mint; 
end

params.rssipos = rssipos;
params.gaussparams = gaussparams;


%{
T = [];
for ki=1:length(keys)
   d = rssidata{ki};
   T = union(T, d(:,1));  
end
T = T / 1000;
%}
T = [min(gtdata_t):params.dt:max(gtdata_t)+params.dt]';

Y = nan*ones(length(keys),length(T));
for ki=1:length(keys)
    d = rssidata{ki};
    if ( size(d,1) > 1 )
        Y(ki,:) = interp1(d(:,1)/1000,d(:,2),T,'nearest');
    end
    %for ti=1:size(d,1)
    %    Y(ki,d(ti,1)/1000>T-params.dt & d(ti,1)/1000<=T) = d(ti,2);
    %end
end


GTinterp = zeros(2,length(T));
GTinterp(1,:) = interp1(gtdata_t,gtdata_x(:,1),T);
GTinterp(2,:) = interp1(gtdata_t,gtdata_x(:,2),T);
% Fake sensor data (simulation)
for ki=1:length(keys)
    for ti=1:size(T)
        gp = params.gaussparams{ki};
        dist = GTinterp(:,ti)-params.rssipos{ki}';
        dist = sqrt(dist(1,:).^2 + dist(2,:).^2);
        dist = dist / 10;
        muint = interp1(gp(:,1),gp(:,2),dist,'linear','extrap');
        sigint = interp1(gp(:,1),gp(:,3),dist,'nearest','extrap');
        Y(ki,ti) = muint+sigint*randn([1 1]);
        %Y(ki,ti) = muint;
    end
end


% Fudge it up
%{
Y(~isnan(Y)) = Y(~isnan(Y)) + 5; % more RSSI
for ki=1:length(keys)
   params.gaussparams{ki}(:,3) = params.gaussparams{ki}(:,3); % more sensor noise
end
%}

%% Uniform Stupid Sampling

nSamples = 10000;
Xest = zeros([2 nSamples]); 


% Create uniform samples

badidx = 1:nSamples;
while (~isempty(badidx))
    fprintf(['Creating initial positions, size = ' num2str(length(badidx)) '\n']);
    Xest(:,badidx) = rand([2 length(badidx)]);
    Xest(1,badidx) = Xest(1,badidx) * (imgbounds(1,2) - imgbounds(1,1)) + imgbounds(1,1);
    Xest(2,badidx) = Xest(2,badidx) * (imgbounds(2,2) - imgbounds(2,1)) + imgbounds(2,1);
    
    [si,sj] = nearestHex(Xest(1,badidx),Xest(2,badidx),params.grid);
    badidx2 = [];
    for i=1:length(si)
        siv = si(i);
        sjv = sj(i);
        if (sjv < 1 || siv < 1 || siv > size(obsmap_hex,2) || sjv > size(obsmap_hex,1) || ...
                obsmap_hex(sjv,siv) > 0.2)
           badidx2 = [badidx2 badidx(i)]; 
        end
    end
    badidx = badidx2;
end


% create fixed samples
%Xest = repmat(gtdata_x(1,:)',[1,nSamples]);

West = ones(nSamples,1) * 1/nSamples;

% Transition function
%transfunc = @(X,dt) X + dt*wSpeed*2*(rand(size(X))-0.5);
% Observation function
%observefunc = @(Y,X) normpdf(Y(1),X(1,:),wifiSigma).*normpdf(Y(2),X(2,:),wifiSigma);
% do the filtering
params.framedir = 'frames/';
params.gtdata_t = gtdata_t;
params.gtdata_x = gtdata_x;

[Xest,West,Xtraj] = SIRFilter(Xest,West,Y,T, @transitionFunc, @observeFunc, params);
%% Visualization
close all;

figure;
subplot(2,1,1); hold on; plot(gtdata_t,gtdata_x(:,1),'+-b'); plot(T,Xtraj(1,:),'o-r'); ylabel('x pos (cm)');
subplot(2,1,2); hold on; plot(gtdata_t,gtdata_x(:,2),'+-b'); plot(T,Xtraj(2,:),'o-r'); ylabel('y pos (cm)');

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
plot(Xtraj(1,:),Xtraj(2,:),'o-r');


