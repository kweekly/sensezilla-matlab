%%
imgbounds = [[-160.11,468.5];[-308.07,159.54]];
%{
[obsmap, entrymap, stopmap] = readFloorplan('floorplan.png');
grid = hexGrid(imgbounds(1,:),imgbounds(2,:),2,2);
obsmap_hex = rect2hexgrid(obsmap, grid);
entrymap_hex = rect2hexgrid(entrymap, grid);
stopmap_hex = rect2hexgrid(stopmap, grid);
%}
%%
dat = importdata('data/test.csv');
data = dat.data(dat.data(:,2)==1,:);
data(:,4) = -data(:,4);
t = data(:,1);
x = data(:,3);
y = data(:,4);

%%
%{
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
hexPlot(visi,visj,'r',grid);
%}
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
[Xest,West,Xtraj] = SIRFilter(Xest,West,[x y]',t, @transitionFunc, @observeFunc);
figure;
subplot(2,1,1); hold on; plot(t,x,'+-b'); plot(t,Xtraj(1,:),'o-g');
subplot(2,1,2); hold on; plot(t,y,'+-b'); plot(t,Xtraj(2,:),'o-g');
figure; hold on;
plot3(Xest(1,:),Xest(2,:),West,'.g'); plot3(x(end),y(end),max(West),'+b');

