%% STEP I Data Loading and Processing
% Copied from EDA_ACme_Zhaoyi.m
clear all;
close all;

% original data has 1sec sampling rate
s1 = importdata('stepI_8d8.csv', ',', 13);
s2 = importdata('fout.csv', ',', 13);
%s2 = importdata('stepI_8e3.csv', ',', 13);
s3 = importdata('stepI_8e4.csv', ',', 13);
s4 = importdata('stepI_937.csv', ',', 13);
DevNum = 2;     % number of devices (assume that we already know)
DevName = {'Overall Load','DELL Monitor','DELL Desktop','IBM Laptop'};

% Some time point is missing, so we gonna fill them with NaN values
TimeF = s1.data(1,1);
TimeL = s1.data(end,1);
T = (TimeL-TimeF)/1000;
RES = nan(T, 4);    % store the full data
KEY = s1.data(:, 1);
RES((KEY-TimeF)/1000+1, 1) = s1.data(:, 2);
KEY = s2.data(:, 1);
RES((KEY-TimeF)/1000+1, 2) = s2.data(:, 2);
KEY = s3.data(:, 1);
RES((KEY-TimeF)/1000+1, 3) = s3.data(:, 2);
KEY = s4.data(:, 1);
RES((KEY-TimeF)/1000+1, 4) = s4.data(:, 2);

% The following codes try to re-sample the data by 10sec rate, and the data
% "RES_samp" stores the result
RES_decimate = 10;
T_samp = ceil(T/RES_decimate);
RES_samp = zeros(T_samp, 4);
for i = 1:DevNum
    LS = 1:RES_decimate:((T_samp-1)*RES_decimate + 1);
    RES_samp(:, i) = RES(LS, i);
end;

%% STEP Ib Min Filter and spike filter
close all;
%{
AGG_filt2 = Kevin_MINFILT(AGG);

plot(AGGX,AGG,AGGX,AGG_filt2);
title(['MIN filtering']);
AGG = AGG_filt2;
%}
figure;
% Plot the sampled results
for i = 1:DevNum
    fprintf(['\nFiltering ' DevName{i} '\n']);
    RES_filt(:,i) = Kevin_MINFILT(RES_samp(:,i));
    subplot(DevNum,1,i);
    plot(RES_samp(:,i),'-b');    
    %plot((1:T_samp),RES_samp(:,i),'-b',(1:T_samp),RES_filt(:,i),'-r');    
    axis tight;	grid on; ylim([0,2.36e5]);
    title(DevName(i));
end;

%% STEP Ic Low-pass filtering
%{
Fwin = 10;
AGG_filt = filter(ones(1,Fwin)/Fwin, 1, AGG, ones(1,Fwin-1)*AGG(1)); 
MSE = sum((AGG_filt - AGG).^2);
fprintf(['\nMSE from filtering = ' num2str(MSE,'%.4e') '\n']);
AGG = AGG_filt;

figure;
hold on;
plot(AGGX,AGG,AGGX,AGG_filt);
title('Low-pass filtering');
%}

%% STEP II Chop-up aggregate data

for i = 1:DevNum
    AGG = RES_filt(:,i)';
    AGGX = (1:T_samp)/60;
    N = length(AGG);


    figure;hold on;
    DERIV_decimate = 20;
    DERIV = [0 diff(AGG(1:DERIV_decimate:N))];
    AGGX_decimated = AGGX(1:DERIV_decimate:N);

    plot(AGGX,AGG,'-b',AGGX_decimated,DERIV,'+r');

    THRESH = 1.3e4;

    Trans = DERIV;
    Trans(abs(Trans)<THRESH) = nan;     
    stem(AGGX_decimated,Trans,'g');
    title(['Chopped up by derivative ' DevName{i}]);
end
%% STEP IIb Iterative Linear Approximation
%{
posA = 1;
THRESH = 1e4;
SEGS = {};
ERRS = {};
while(posA < length(AGG))
    posB = posA+1;
    err_accum = [0];
    while(posB < length(AGG))
       X = posA:posB;
       m = (AGG(posB)-AGG(posA))/(posB-posA);
       Y = m*X + AGG(posA) - m*posA;
       err = sum(Y - AGG(X)) / (posB-posA);
       err_accum(end+1) = err;
       
       %{
       Y = ones(1,posB-posA+1)*mean(AGG(X));
       err = sum(abs(Y - AGG(X))) / (posB-posA);
       err_accum(end+1) = err;
       %}
       %{
       close all;
       figure; hold on;   
       plot(Y,'r');
       plot(AGG(X)     ,'b');
       plot(cumsum(abs(Y - AGG(X))),'g');
       pause;
       %}
       
       if (abs(err) > THRESH)
           SEGS{end+1} = [posA posB];
           ERRS{end+1} = err_accum;
           posA = posB;
           break;
       end
       posB = posB + 1;
    end
    if ( posB == length(AGG) )
       SEGS{end+1} = [posA posB];
       err_accum(end+1) = err_accum(end);
       ERRS{end+1} = err_accum;
       posA = posB;
       break;
    end
end

figure; hold on;
plot(AGGX,AGG,'-b');
for i=1:length(SEGS)
   AVG = mean(AGG(SEGS{i}(1):SEGS{i}(2)));
   plot(AGGX(SEGS{i}),AGG(SEGS{i}),'Color','r','LineWidth',2,'Marker','o');
   %plot(AGGX(SEGS{i}),ones(1,2)*AVG,'Color','r','LineWidth',2,'Marker','o');
   
   plot(AGGX(SEGS{i}(1):SEGS{i}(2)),ERRS{i},'Color','g','LineWidth',2);
end
%}