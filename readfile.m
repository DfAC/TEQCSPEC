function [t_samp mjl SAT n] = readfile(N,n,i,A,file)
% See also TEQCSPEC, CHECKFILE, SCANFILE
%
% History
% 22 Feb 2009 created using Matlab R2008b

global Sat_Capacity

SAT(1:N,1:Sat_Capacity)=NaN;  % 1-32 : GPS ; 33-64 : GLONASS
k=0;
fid=fopen(file);

sats=A;

% h=waitbar(0,['Please wait... Importing ' char(file)]);

while 1;
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    k=k+1;
    if k==i-2, t_samp=char(tline);end %% sample rate
    if k==i-1, mjl=char(tline);end    %% MJL
    if k > i,
%         waitbar(k/N,h);drawnow
        if mod(k,5000)==0
            disp(strcat(num2str(100*k/N),'%'))
        end
        n=n+1;
        if sats==0
            A=UpdateSatInview(tline);
        else
            A=str2num(tline);
        end
        if sats~=0;
            SAT(n,sats)=A;
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            k=k+1;
            A=UpdateSatInview(tline);
        end
        if A~=-1
            sats=A;
        end
    end
end
fclose(fid);
% close(h);

