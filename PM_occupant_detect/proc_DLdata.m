%% Loading data
close all;
datdir = ['D:\\documents\\ucb\\singapore\\data\\PMsensortest\\SWARM_02282013'];
%datdir = ['D:\\documents\\ucb\\singapore\\data\\PMsensortest\\SWARM_03072013'];
pmdatafile = [datdir '\\data.mat'];
camdatafile = [datdir '\\tfound.mat'];
load(pmdatafile);

fprintf(['Loaded PM data from ' pmdatafile '\n']);

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
load(camdatafile);
fprintf(['Loaded camera data from ' camdatafile '\n']);
%% Window filtering

twind = 300;
ssratio = 30;
ymax = 4e-3; 

fprintf(['Parameters: twind=' num2str(twind) ' ssratio=' num2str(ssratio) ' Particle Size=2um\n']);

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
combf = datf;

subaxis(3);
dat = ch1_2u;
doaplot_DLdata; 
r = corrcoef(tfoundf,datf);
title(['Ch1 >2um particles. ' num2str(twind) 's Time window r=' num2str(r(1,2))]); ylim([0 ymax]);
ch1f = datf;

subaxis(4);
dat = ch2_2u;
doaplot_DLdata; 
r = corrcoef(tfoundf,datf);
title(['Ch2>2um particles. ' num2str(twind) 's Time window r=' num2str(r(1,2))]);  ylim([0 ymax]);
ch2f = datf;

subaxis(5);
dat = ch3_2u;
doaplot_DLdata; 
r = corrcoef(tfoundf,datf);
title(['Ch3 >2um particles. ' num2str(twind) 's Time window r=' num2str(r(1,2))]);  ylim([0 ymax]);
ch3f = datf;

subaxis(6);
dat = ch4_2u;
doaplot_DLdata; 
r = corrcoef(tfoundf,datf);
title(['Ch4 >2um particles.' num2str(twind) 's Time window r=' num2str(r(1,2))]);  ylim([0 ymax]);
ch4f = datf;

subaxis(7);
dat = ch5_2u;
doaplot_DLdata; 
r = corrcoef(tfoundf,datf);
title(['Ch5 >2um particles.' num2str(twind) 's Time window r=' num2str(r(1,2))]);  ylim([0 ymax]);
ch5f = datf;

%% Perform calibration
% for DSM501A 2um output
cal_a = 8.475e4;
cal_b = 45.71;

combf = combf * cal_a + cal_b; % in units of pcs/L

%%  Finding shift amount 

[i,j1] = max(xcorr( tfoundf, combf ));
[i,j2] = max(xcorr( tfoundf, tfoundf));
shift = j2 - j1;
tfoundf_sh = tfoundf(1:end-shift);
combf_sh = combf(1+shift:end);
tf_sh = tf(1+shift:end);
rold = corrcoef(tfoundf,combf);
rnew = corrcoef(tfoundf_sh,combf_sh);
fprintf(['Detected shift amount=' num2str(shift) ' (' num2str(shift * (twind/ssratio)) 's) old r=' num2str(rold(1,2)) ' new r=' num2str(rnew(1,2)) '\n']);

%% Numerical Probability Distributions
xbins = [0:4e-3:2.8e-2 inf]; % x are in units of "camera occurances per second"
ybins = [40:20:200 inf]; % y are in units of pcs/L
figure;
[Z, counts, xbins, ybins] = histheatmap(tfoundf_sh,combf_sh,xbins,ybins);
%yscale = 1e-3;
yscale = 1;
fprintf('\nObserved Occurances Pr(x,y):\n            ');
for j=1:length(ybins)-1
   fprintf('& %3.1f - %3.1f ',ybins(j)/yscale,ybins(j+1)/yscale);
end
fprintf('\\\\\\hline\\hline\n');

xscale = 1e-2;
for i=1:length(xbins)-1
   fprintf('%4.1f - %4.1f',xbins(i)/xscale,xbins(i+1)/xscale);
   for j=1:length(ybins)-1
      fprintf(' & %9d%',Z(j,i));
   end
   fprintf(' \\\\\\hhline{~*{11}{-}}\n');
end

%% Calculate Binary Detector

PM_th = 100;
Cam_th = 1e-2;

TP= nnz(tfoundf_sh < Cam_th & combf_sh < PM_th)/nnz(combf_sh < PM_th);
FP= nnz(tfoundf_sh >= Cam_th & combf_sh < PM_th)/nnz(combf_sh < PM_th);
TN= nnz(tfoundf_sh >= Cam_th & combf_sh >= PM_th)/nnz(combf_sh >= PM_th);
FN= nnz(tfoundf_sh < Cam_th & combf_sh >= PM_th)/nnz(combf_sh >= PM_th);
TOTN = nnz((tfoundf_sh >= Cam_th & combf_sh < PM_th) | (tfoundf_sh < Cam_th & combf_sh >= PM_th) )/length(combf_sh);

fprintf(['TRUE  POS: Pr(Cam < ' num2str(Cam_th) ' | PM_th < ' num2str(PM_th) ')=' num2str( TP ) '\n']);
fprintf(['FALSE POS: Pr(Cam >= ' num2str(Cam_th) ' | PM_th < ' num2str(PM_th) ')=' num2str( FP ) '\n']);

fprintf(['TRUE  NEG: Pr(Cam >= ' num2str(Cam_th) ' | PM_th >= ' num2str(PM_th) ')=' num2str( TN ) '\n']);
fprintf(['FALSE NEG: Pr(Cam < ' num2str(Cam_th) ' | PM_th >= ' num2str(PM_th) ')=' num2str(  FN ) '\n']);

fprintf(['TOTAL NEG: ' num2str( TOTN ) '\n']);

fprintf('\n');

%% Calculate Pearson's Coeff. for coarseness of particles
f = figure;
dat = mean([ch6_05u ch7_05u ch8_05u],2);
doaplot_DLdata;
r05u = corrcoef(tfoundf(1:end-shift),datf(1+shift:end));

dat = mean([ch1_1u ch2_1u ch3_1u ch4_1u ch5_1u],2);
doaplot_DLdata;
r1u = corrcoef(tfoundf(1:end-shift),datf(1+shift:end));

dat = mean([ch1_2u ch2_2u ch3_2u ch4_2u ch5_2u],2);
doaplot_DLdata;
r2u = corrcoef(tfoundf(1:end-shift),datf(1+shift:end));

rcam = corrcoef(tfoundf,tfoundf);
fprintf('Correlation Coefficients- 0.5u:%.4f 1u:%.4f 2u:%.4f cam:%.4f\n',r05u(1,2),r1u(1,2),r2u(1,2),rcam(1,2));
close(f);

%% Kalman Filter
pfitrev = polyfit(combf_sh, tfoundf_sh, 1);
pfit = polyfit(tfoundf_sh, combf_sh, 1);
pfit_est = pfitrev(1)*combf_sh + pfitrev(2);

% From: http://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf
% transition matrix 
decayO  = 0.001;
A = eye(2) + [[-decayO 0];[0 0]];
% covariance of model noise
maxdiffO = 10 * max(combf_sh);
maxdiffD = 0.1 * mean(combf_sh);
Q = [[maxdiffO^2 0];[0 maxdiffD^2]];

% output matrix
H = [pfit(1) 1];
% covariance of sensor noise
R = 0.001 * max(combf_sh);

x_est = zeros([2 length(combf_sh)]);
x_est(:,1) = [combf_sh(1)/2;combf_sh(1)/2];
P = 100 * ones([2 2]);
for i=2:length(combf_sh)
    x_est(:,i) = A*x_est(:,i-1);
    P = A*P*A' + Q;
    
    K = P*H'*(H*P*H'+R)^(-1);
    x_est(:,i) = x_est(:,i) + K*(combf_sh(i) - H*x_est(:,i));
    P = (eye(2) - K*H)*P;
end

fprintf('LSE linfit=%.3f  kalman=%.3f\n',norm((pfitrev(1)*combf_sh+pfitrev(2)) - tfoundf_sh),norm(x_est(1,:)' - tfoundf_sh));

figure;
subaxis(3,1,1,  'Spacing', 0.03, 'Padding', .05, 'MarginRight',0.01,'MarginLeft',0.08,'MarginBottom',0.1,'MarginTop',0.05); 
hold on;
plot(tf_sh,tfoundf_sh,'m'); 
plot(tf_sh,pfit_est,'g');
plot(tf_sh,x_est(1,:)); title('Estimate of Local Occupancy');
subaxis(2);
plot(tf_sh,x_est(2,:)); title('Estimate of Room Occupancy'); 
subaxis(3); hold on;
plot(tf_sh,H*x_est);
plot(tf_sh,combf_sh,'r'); title('Comparison est. output vs. real');



%% "For Show" Plots

figure; 
subaxis(2,1,1,  'Spacing', 0.03, 'Padding', 0, 'MarginRight',0.01,'MarginLeft',0.08,'MarginBottom',0.1,'MarginTop',0.05); 
hold on;
stem(tfound/3600,max(tfoundf_sh/xscale)*ones([1 length(tfound)]),'g.');
plot(tf_sh/3600,tfoundf_sh/xscale,'m'); ylabel({'Camera Occurances';'per 100s'});
set(gca,'XTickLabel','');
set(gca,'XTick',[]);
axis tight;
subaxis(2);
plot(tf_sh/3600,combf_sh/yscale); xlabel('Hours'); ylabel({'Concentration of';'\geq 2 \mum particles (per L)'});
axis tight;

%{
figure; hold on;
plot(tfoundf_sh/xscale,  combf_sh/yscale, '.'); xlabel({'Camera Occurances';'per 100s'}); ylabel({'PM Concentration';'(2 \mum particles/L)'});
tvar = [0:1e-3:max(tfoundf_sh)];
plot(tvar/xscale,(tvar * pfit(1) + pfit(2))/yscale,'Color','m','LineWidth',2);
axis tight;
%}

figure; hold on;
plot(combf_sh/yscale,  tfoundf_sh/xscale, '.'); ylabel({'Camera Occurances';'per 100s'}); xlabel({'Concentration of';'\geq 2 \mum particles (per L)'});
tvar = [0 max(combf_sh/yscale)];
plot(tvar/yscale,(tvar * pfitrev(1) + pfitrev(2))/xscale,'Color','m','LineWidth',2);
axis tight;


