function [tfound,vidobj] = proc_landmarks(fname, lmarks)

%%
vidobj = mmreader(fname);
h = vidobj.Height;
w = vidobj.Width;
frameskip = 5;
frames = vidobj.NumberOfFrames / frameskip;

lmark_xblk = -10:10;
lmark_yblk = -20:20;

lmark_hue = 0.01459;
hue_th = 0.02;
lmark_sat = 1;
sat_th = 0.1;
lmark_val = 0.7;
val_th = 0.4;

%%
lmres = zeros(frames,size(lmarks,1));
for frameno=1:frames
   imgmat = read(vidobj, frameno * frameskip);
  
   for lmi=1:size(lmarks,1)
        imghsv = rgb2hsv(imgmat(lmarks(lmi,2)+lmark_yblk,lmarks(lmi,1)+lmark_xblk,:));
        hdiff = imghsv(:,:,1) - lmark_hue;
        img_th = (abs(mod(hdiff + 0.5,1)-0.5) <= hue_th) & ...
            (abs(imghsv(:,:,2) - lmark_sat) <= sat_th) & ...
            (abs(imghsv(:,:,3) - lmark_val) <= val_th);
        lmres(frameno,lmi) = sum(sum(img_th));
   end
   if mod(frameno,round(vidobj.FrameRate)*10) == 0
       fprintf(['Processing ' num2str(round(frameno/vidobj.FrameRate * frameskip)) ' of ' num2str(round(frames/vidobj.FrameRate * frameskip)) ' seconds\n']);
   end
end


%%
lmres_th = lmres > (length(lmark_xblk)*length(lmark_yblk)/2);

t = (0:length(lmres)-1) / vidobj.FrameRate * frameskip;
fidx = find(sum(lmres_th,2) < size(lmarks,1) - 2); % Find at least 2 landmarks blocked
if length(fidx) > 0
    tfound = t(fidx(1));
    lastt = tfound;
    tblockn = 1;
    tblock_th = 0.5;
    for fii=1:length(fidx)
       ti = fidx(fii);
       if t(ti) - lastt < tblock_th
          tfound(end) = tfound(end) + t(ti);
          tblockn = tblockn + 1;
          lastt = t(ti);
       else
           tfound(end) = tfound(end) / tblockn;
           tfound = [tfound;t(ti)];
           tblockn = 1;
           lastt = t(ti);
       end
    end

    tfound(end) = tfound(end) / tblockn;
else
    tfound = zeros(0,1);
end

%{
idx = 0;
last_row = ones(1,size(lmarks,1));
while idx < frames
    % find changes
   idx = idx + find( lmres_th(idx + 1:end,:) ~= repmat(last_row,size(lmres_th,1) - idx,1) , 1 );
   row = lmres_th(idx,:);
   % check which direction things are moving
   if sum(last_row(1:2)) < 2 && sum(last_row(end-1:end)) == 2 && sum(row(end-1:end)) < 2
      fprintf(['To the right now yall t=' num2str(t(idx)) '\n']);
   elseif sum(last_row(end-1:end)) < 2 && sum(last_row(1:2)) == 2 && sum(row(1:2)) < 2
      fprintf(['To the left now yall t=' num2str(t(idx)) '\n']); 
   end
   last_row = row;
end
%}
end 