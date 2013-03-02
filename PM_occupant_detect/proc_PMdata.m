%%
dat = importdata('CREST_02262013\pmdata.par');
t = dat(:,1)/1000.0;
v = dat(:,2);
while ( v(1) == 1 || (v(1) == 0 && v(2) == 0) )
   t = t(2:end);
   v = v(2:end);
end

PW = t(2:2:end) - t(1:2:end-1);
PWs = t(1:2:end-1);
PWe = t(2:2:end);

twind = 1;
tres = 1;
T = t(1) + twind/2:tres:t(end)-twind/2;
PM = zeros(size(T));
PWA = zeros(size(T));
PC = zeros(size(T));
for TI=1:length(T)
   t1 = T(TI) - twind/2;
   t2 = T(TI) + twind/2;
   PM(TI) = sum(PW .* ((PWs >= t1) & (PWe <= t2)) + (PWe - t1) .* ((PWs <= t1) & (PWe > t1)) + (t2-PWs) .* ((PWs < t2) & (PWe >= t2)) ) / twind;
   PC(TI) = sum((PWs >= t1) & (PWs <= t2) | (PWe >= t1) & (PWe <= t2));
   PWA(TI) = PM(TI) * twind / PC(TI);
end


filtamt = 1;
%%

PWAf = filter(ones([1 filtamt])/filtamt,1,PWA);
PMf = filter(ones([1 filtamt])/filtamt,1,PM);
PCf = filter(ones([1 filtamt])/filtamt,1,PC);

figure;
hold on;
stem(tfound,ones(size(tfound)),'LineWidth',2,'Color','m');
plot(T,PWAf,'r');
plot(T,PMf+0.01,'b');
plot(T,PCf/max(PCf) * 0.02,'g');
legend('Camera','Pulse Width','Duty Cycle','Num Pulses');
axis tight; ylim([0,max(PM)]);

look = 600;
PMchunks = [];
%figure; hold on;
for tfi=1:length(tfound)
   tfv = tfound(tfi);
   idx = find(T >= tfv & T < tfv+look);
   Tfl = T(idx) - tfv;
%   PWAfl = PWAf(idx);
   PMfl = PMf(idx);
%   PCfl = PCf(idx);
%   plot(Tfl,PWAfl,'r');
  plot(Tfl,PMfl+0.01,'b');
%   plot(Tfl,PCfl/max(PCf) * 0.02,'g');
   if (size(PMfl,2) < size(PMchunks,2))
      PMfl = [PMf zeros(1,size(PMchunks,2) - size(PMfl,2))];
   else
      PMchunks = [PMchunks zeros(size(PMchunks,1),size(PMfl,2) - size(PMchunks,2))];
   end
   PMchunks = [PMchunks;PMfl];
end


%%
figure;
b = [min(min(PMchunks)) max(max(PMchunks))];
bins = 0:0.15/200:0.15;
[X,Y] = meshgrid(Tfl,bins);
histZ = zeros(length(bins),length(Tfl));

%{
combidx = combnk(1:size(PMchunks,1),15);
combv = zeros(size(combidx,1),size(PMchunks,2));

hold on;
for i=1:size(combidx,1)
    %plot(mean(PMchunks(combidx(i,:),:)));
    combv(i,:) = mean(PMchunks(combidx(i,:),:));
end


for i=1:size(combv,2)
    histZ(:,i) = histc(combv(:,i),bins);
end

%}

for i=1:size(PMchunks,2)
   histZ(:,i) = histc(bootstrp(1000,@mean,PMchunks(:,i)) , bins);
end

surf(X,Y,histZ,'EdgeColor','none'); ylim([0 0.15]); view([0 90]);


