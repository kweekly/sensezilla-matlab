[obsmap, entrymap, stopmap] = readFloorplan('floorplan.png');
imgbounds = [[-160.11,468.5];[-308.07,159.54]];
grid = hexGrid(imgbounds(1,:),imgbounds(2,:),1);
obsmap_hex = rect2hexgrid(obsmap, grid);
entrymap_hex = rect2hexgrid(entrymap, grid);
stopmap_hex = rect2hexgrid(stopmap, grid);

figure;
hold on;
surf(grid.xs,grid.ys,obsmap_hex-5,'EdgeColor','None'); view([0, 90]);
dat = importdata('data/test.csv');
data = dat.data(dat.data(:,2)==1,:);
plot(data(:,3),-data(:,4),'g-+');
