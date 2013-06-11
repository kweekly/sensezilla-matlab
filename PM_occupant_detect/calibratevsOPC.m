%% Loading data
close all;
%datdir = ['D:\\documents\\ucb\\singapore\\data\\PMsensortest\\SWARM_02282013'];
datdir = ['D:\\documents\\ucb\\singapore\\data\\PMsensortest\\SWARM_03072013'];
pmdatafile = [datdir '\\data.mat'];
opcdatafile = [datdir '\\OPCdata.mat'];

load(pmdatafile);
% Every sample is approximately 50ms

pmdata.t = data(:,1) / 2000;
pmdata.t = pmdata.t - pmdata.t(1);

% turn into "low times"
data(:,3:15) = 1-data(:,3:15)./repmat(data(:,2),1,15-3+1);

pmdata.ch1_2u = data(:,3);
pmdata.ch1_1u = data(:,4);
pmdata.ch2_2u = data(:,5);
pmdata.ch2_1u = data(:,6);
pmdata.ch3_2u = data(:,7);
pmdata.ch3_1u = data(:,8);
pmdata.ch4_2u = data(:,9);
pmdata.ch4_1u = data(:,10);
pmdata.ch5_2u = data(:,11);
pmdata.ch5_1u = data(:,12);

pmdata.ch6_05u = data(:,13);
pmdata.ch7_05u = data(:,14);
pmdata.ch8_05u = data(:,15);

pmdata.mean_2u = mean([pmdata.ch1_2u pmdata.ch2_2u pmdata.ch3_2u pmdata.ch4_2u pmdata.ch5_2u],2);
pmdata.mean_1u = mean([pmdata.ch1_1u pmdata.ch2_1u pmdata.ch3_1u pmdata.ch4_1u pmdata.ch5_1u],2);
pmdata.mean_05u = mean([pmdata.ch6_05u pmdata.ch7_05u pmdata.ch8_05u],2);

clear data;

opcdata = load(opcdatafile);
opcdata.t = opcdata.t - opcdata.t(1);

%% Averaging

fprintf('ch1 2u averaging \n');
pmdata.ch1_2u_cal = do_calibrate_sametime(pmdata.t, pmdata.ch1_2u, opcdata.t, opcdata.pm2);
fprintf('ch2 2u averaging \n');
pmdata.ch2_2u_cal = do_calibrate_sametime(pmdata.t, pmdata.ch2_2u, opcdata.t, opcdata.pm2);
fprintf('ch3 2u averaging \n');
pmdata.ch3_2u_cal = do_calibrate_sametime(pmdata.t, pmdata.ch3_2u, opcdata.t, opcdata.pm2);
fprintf('ch4 2u averaging \n');
pmdata.ch4_2u_cal = do_calibrate_sametime(pmdata.t, pmdata.ch4_2u, opcdata.t, opcdata.pm2);
fprintf('ch5 2u averaging \n');
pmdata.ch5_2u_cal = do_calibrate_sametime(pmdata.t, pmdata.ch5_2u, opcdata.t, opcdata.pm2);
fprintf('mean 2u averaging \n');
pmdata.mean_2u_cal = do_calibrate_sametime(pmdata.t, pmdata.mean_2u, opcdata.t, opcdata.pm2);

fprintf('ch1 1u averaging \n');
pmdata.ch1_1u_cal = do_calibrate_sametime(pmdata.t, pmdata.ch1_1u, opcdata.t, opcdata.pm1);
fprintf('ch2 1u averaging \n');
pmdata.ch2_1u_cal = do_calibrate_sametime(pmdata.t, pmdata.ch2_1u, opcdata.t, opcdata.pm1);
fprintf('ch3 1u averaging \n');
pmdata.ch3_1u_cal = do_calibrate_sametime(pmdata.t, pmdata.ch3_1u, opcdata.t, opcdata.pm1);
fprintf('ch4 1u averaging \n');
pmdata.ch4_1u_cal = do_calibrate_sametime(pmdata.t, pmdata.ch4_1u, opcdata.t, opcdata.pm1);
fprintf('ch5 1u averaging \n');
pmdata.ch5_1u_cal = do_calibrate_sametime(pmdata.t, pmdata.ch5_1u, opcdata.t, opcdata.pm1);
fprintf('mean 1u averaging \n');
pmdata.mean_1u_cal = do_calibrate_sametime(pmdata.t, pmdata.mean_1u, opcdata.t, opcdata.pm1);

fprintf('ch6 05u averaging \n');
pmdata.ch6_05u_cal = do_calibrate_sametime(pmdata.t, pmdata.ch6_05u, opcdata.t, opcdata.pm05);
fprintf('ch7 05u averaging \n');
pmdata.ch7_05u_cal = do_calibrate_sametime(pmdata.t, pmdata.ch7_05u, opcdata.t, opcdata.pm05);
fprintf('ch8 05u averaging \n');
pmdata.ch8_05u_cal = do_calibrate_sametime(pmdata.t, pmdata.ch8_05u, opcdata.t, opcdata.pm05);
fprintf('mean 05u averaging \n');
pmdata.mean_05u_cal = do_calibrate_sametime(pmdata.t, pmdata.mean_05u, opcdata.t, opcdata.pm05);

%% Sliding Window filtering to 5m
nsamp = 10;

pmdata.ch1_2u_5m = do_calibrate_filter(pmdata.ch1_2u_cal,nsamp);
pmdata.ch2_2u_5m = do_calibrate_filter(pmdata.ch2_2u_cal,nsamp);
pmdata.ch3_2u_5m = do_calibrate_filter(pmdata.ch3_2u_cal,nsamp);
pmdata.ch4_2u_5m = do_calibrate_filter(pmdata.ch4_2u_cal,nsamp);
pmdata.ch5_2u_5m = do_calibrate_filter(pmdata.ch5_2u_cal,nsamp);
pmdata.mean_2u_5m = do_calibrate_filter(pmdata.mean_2u_cal,nsamp);

pmdata.ch1_1u_5m = do_calibrate_filter(pmdata.ch1_1u_cal,nsamp);
pmdata.ch2_1u_5m = do_calibrate_filter(pmdata.ch2_1u_cal,nsamp);
pmdata.ch3_1u_5m = do_calibrate_filter(pmdata.ch3_1u_cal,nsamp);
pmdata.ch4_1u_5m = do_calibrate_filter(pmdata.ch4_1u_cal,nsamp);
pmdata.ch5_1u_5m = do_calibrate_filter(pmdata.ch5_1u_cal,nsamp);
pmdata.mean_1u_5m = do_calibrate_filter(pmdata.mean_1u_cal,nsamp);

pmdata.ch6_05u_5m = do_calibrate_filter(pmdata.ch6_05u_cal,nsamp);
pmdata.ch7_05u_5m = do_calibrate_filter(pmdata.ch7_05u_cal,nsamp);
pmdata.ch8_05u_5m = do_calibrate_filter(pmdata.ch8_05u_cal,nsamp);
pmdata.mean_05u_5m = do_calibrate_filter(pmdata.mean_05u_cal,nsamp);

opcdata.pm03_5m = do_calibrate_filter(opcdata.pm03,nsamp);
opcdata.pm05_5m = do_calibrate_filter(opcdata.pm05,nsamp);
opcdata.pm07_5m = do_calibrate_filter(opcdata.pm07,nsamp);
opcdata.pm1_5m = do_calibrate_filter(opcdata.pm1,nsamp);
opcdata.pm2_5m = do_calibrate_filter(opcdata.pm2,nsamp);
opcdata.pm5_5m = do_calibrate_filter(opcdata.pm5,nsamp);


%% Remove some nastiness
idx1 = 1356;
idx2 = 2044;

opcdata.t_rn = [opcdata.t(1:idx1);opcdata.t(idx2:end)];

opcdata.pm2_rn = [opcdata.pm2_5m(1:idx1);opcdata.pm2_5m(idx2:end)];
opcdata.pm1_rn = [opcdata.pm1_5m(1:idx1);opcdata.pm1_5m(idx2:end)];
opcdata.pm05_rn = [opcdata.pm05_5m(1:idx1);opcdata.pm05_5m(idx2:end)];

pmdata.mean_2u_rn = [pmdata.mean_2u_5m(1:idx1);pmdata.mean_2u_5m(idx2:end)];
pmdata.mean_1u_rn = [pmdata.mean_1u_5m(1:idx1);pmdata.mean_1u_5m(idx2:end)];
pmdata.mean_05u_rn = [pmdata.mean_05u_5m(1:idx1);pmdata.mean_05u_5m(idx2:end)];

%% Pretty plots
% Calculate some values
coeff2u = polyfit(opcdata.pm2_5m, pmdata.mean_2u_5m, 1);
coeff2u_rev = polyfit(pmdata.mean_2u_5m, opcdata.pm2_5m, 1);
coeff1u = polyfit(opcdata.pm1_5m, pmdata.mean_1u_5m, 1);
coeff1u_rev = polyfit(pmdata.mean_1u_5m, opcdata.pm1_5m, 1);
coeff05u = polyfit(opcdata.pm05_5m, pmdata.mean_05u_5m, 1);
coeff05u_rev = polyfit(pmdata.mean_05u_5m, opcdata.pm05_5m, 1);

mse2 = sqrt( mean( (pmdata.mean_2u_5m * coeff2u_rev(1) + coeff2u_rev(2) - opcdata.pm2_5m).^2 ) ) / ( max(opcdata.pm2_5m) - min(opcdata.pm2_5m))
mse1 = sqrt( mean( (pmdata.mean_1u_5m * coeff1u_rev(1) + coeff1u_rev(2) - opcdata.pm1_5m).^2) )  / ( max(opcdata.pm1_5m) - min(opcdata.pm1_5m))
mse05 = sqrt( mean( (pmdata.mean_05u_5m * coeff05u_rev(1) + coeff05u_rev(2) - opcdata.pm05_5m).^2) )  / ( max(opcdata.pm05_5m) - min(opcdata.pm05_5m))


%% Scatter 2u
figure; hold on;
plot(pmdata.mean_2u_5m*1e3,opcdata.pm2_5m,'.');
plot([0 max(pmdata.mean_2u_5m*1e3)],(coeff2u_rev(2) + [0 coeff2u_rev(1)*max(pmdata.mean_2u_5m)]),'LineWidth',2,'Color','m');
%ylim([0 2.5]);
axis tight;
ylim([0 250]);
xlabel('Mean PM Sensor Ratio (parts per 1000)');
ylabel('Concentration of \geq 2\mum particles (pcs/liter)');


