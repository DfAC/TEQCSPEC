function L = checkfile(file)
fid=fopen(file);
L=0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    L = L + 1;
end
fclose(fid);

%fid=fopen('Ball1830.mp1');
%ln=0;
%while 1
%    ln = ln+1;
%    tline = fgetl(fid);
%    if ~ischar(tline), break, end
%    if ln==5, % get number of satellites in first header
%        disp(tline); 
%        sats=str2num(tline);
%    end
%end
%fclose(fid);
