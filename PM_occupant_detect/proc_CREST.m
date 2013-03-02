% for TEST001.AVI
%{
lmarks = [[29,300];
          [134,146];
          [237,296];
          [343,143];
          [436,314]];
%}    

% others For 02/25/2013
lmarks = [[135,299];
          [242,144];
          [342,296];
          [443,156];
          [538,326]];
      
% For 02/26/2013
lmarks = [[156,310];
  [253,157];
  [357,309];
  [460,158];
  [556,323]];

% For SWARM_02282013
lmarks = [[110,379];
            [199,258];
            [276,386];
            [367,268];
            [435,403]];
      
dd = 'SWARM_02282013';
files = dir([dd '/*.AVI']);
toff = 0;
tfound = zeros(0,1);
for fi=1:length(files)
    fname = [dd '/' files(fi).name];
   fprintf([fname '\n']);
   if fi >= 2
      toff = toff + lasto.Duration 
   end
   [tf,lasto] = proc_landmarks(fname,lmarks);
   tfound = [tfound;tf + toff];
end