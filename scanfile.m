function A = scanfile(file,i)
fid=fopen(file);
ln=0;
while 1
    ln=ln+1;
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if ln==i, % get number of satellites in first header
        %disp(tline); 
        A=UpdateSatInview(tline);break;
    end
end
fclose(fid);
