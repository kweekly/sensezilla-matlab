
%for twind=[5 15 30 60 600]
   twindsamp = floor(twind / 0.05 / ssratio);
   datf = zeros([floor(length(dat)/twindsamp) 1]);
   tf = zeros(size(datf));
   tfoundf = zeros(size(tf));
   for i=1:length(datf)
       t1 = t(1+(i-1)*twindsamp);
       t2 = t(1+(i)*twindsamp);
      datf(i) = mean(dat(1+(i-1)*twindsamp:1+i*twindsamp));
      tf(i) = (t1 + t2)/2;
      tfoundf(i) = sum(tfound < t2 & tfound >= t1) / (t2 - t1);
   end
   
   wts = [1/(ssratio*2);repmat(1/ssratio,ssratio-1,1);1/(ssratio*2)];
   datf = conv(datf,wts,'same');
   tfoundf = conv(tfoundf,wts,'same');

   
   %fig = figure; 
   hold on;
   plot(tf,datf); axis tight;
   %stem(tfound,mean(datf)*ones([1 length(tfound)]),'om'); 
   %plot(tf,max(datf)/max(tfoundf) * tfoundf,'m');
   title(['twind=' num2str(twind)]);
%end