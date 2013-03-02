
load('SWARM_02282013\data.mat');

% Every sample is approximately 50ms

t = data(:,1) / 2000;

% turn into "low times"
data(:,3:15) = 1-data(:,3:15)./repmat(data(:,2),1,15-3+1);

ch1_2u = data(:,3);
ch1_1u = data(:,4);
ch2_2u = data(:,5);
ch2_1u = data(:,6);
ch3_2u = data(:,7);
ch3_1u = data(:,8);
ch4_2u = data(:,9);
ch4_1u = data(:,10);
ch5_2u = data(:,11);
ch5_1u = data(:,12);

ch6_05u = data(:,13);
ch7_05u = data(:,14);
ch8_05u = data(:,15);

dat = ch1_2u;
twind = 60;
ssratio = 10;
load('SWARM_02282013\tfound.mat');
%%

twind = 300;
ssratio = 30;
ymax = 4e-3; 

%dat = ch1_1u;
%doaplot_DLdata; title('Ch1 >1um particles. 1m Time window');
f = figure;
doaplot_DLdata;

subaxis(7,1,1,  'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
hold on;
stem(tfound,max(tfoundf)*ones([1 length(tfound)]),'g.');
plot(tf,tfoundf,'m');  xlim([min(tf) max(tf)]);
title('Camera Activity. 1m Time Window');

subaxis(2);
dat = mean([ch1_2u ch2_2u ch3_2u ch4_2u ch5_2u],2);
doaplot_DLdata; 
r = corrcoef(tfoundf,datf);
title(['Comb. >2um particles. ' num2str(twind) 's Time window r=' num2str(r(1,2))]);  ylim([0 ymax]);

subaxis(3);
dat = ch1_2u;
doaplot_DLdata; 
r = corrcoef(tfoundf,datf);
title(['Ch1 >2um particles. ' num2str(twind) 's Time window r=' num2str(r(1,2))]); ylim([0 ymax]);

subaxis(4);
dat = ch2_2u;
doaplot_DLdata; 
r = corrcoef(tfoundf,datf);
title(['Ch2>2um particles. ' num2str(twind) 's Time window r=' num2str(r(1,2))]);  ylim([0 ymax]);

subaxis(5);
dat = ch3_2u;
doaplot_DLdata; 
r = corrcoef(tfoundf,datf);
title(['Ch3 >2um particles. ' num2str(twind) 's Time window r=' num2str(r(1,2))]);  ylim([0 ymax]);

subaxis(6);
dat = ch4_2u;
doaplot_DLdata; 
r = corrcoef(tfoundf,datf);
title(['Ch4 >2um particles.' num2str(twind) 's Time window r=' num2str(r(1,2))]);  ylim([0 ymax]);

subaxis(7);
dat = ch5_2u;
doaplot_DLdata; 
r = corrcoef(tfoundf,datf);
title(['Ch5 >2um particles.' num2str(twind) 's Time window r=' num2str(r(1,2))]);  ylim([0 ymax]);


