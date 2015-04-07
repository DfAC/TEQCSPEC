function out=teqcspec(files);
 close(findall(0,'Type','figure')); 
 clear all
dbstop if error
dbstop at 146
  
%LKB stuff
Constant.PlotName ='QC_'; %see below for name change!
Constant.PlotResolution=300;
  

%TEQCSPEC  TEQC multipath spectrum analysis
%++++++++++++++++++++++++++++++++
%   TEQCSPEC(FILENAME) makes a geographically oriented plot of 
%   multipath in the band 0.05 - 0.06 Hz (default) in a defined
%   time range (see example below).
%
%   Clement.Ogaja@gmail.com
% edit after Lei by LKB (2014-2015)


frqlim = 0.2;  %Max. plotting frequency (Hz)
cutff1 = 0.005; %Lower cut off frequency (Hz)
cutff2 = 0.008; %Upper cut off frequency (Hz)%


% t1=1000;t2=1500; %t1=5600;t2=6380; % for FFT spectra [seconds from 00:00:00 UTC]
%T1=1;T2=6200;%3200;T2=6200; % for polar plots [seconds from 00:00:00 UTC]
% T1=5400;T2=9305;%27204;%7505; % for polar plots [seconds from 00:00:00 UTC]
%t1=4199;t2=5199; % for FFT spectra [seconds from 00:00:00 UTC]
%T1=3199;T2=6199; % for polar plots [seconds from 00:00:00 UTC]
%
%T1=T1+17830;
%T2=T2+17830;
% T1=1;T2=6200;
% t1=t1+17431;
% t2=t2+17431;

% cutff1=0.04;
% cutff2=0.05;

% [t1 t2]=FindBestPeriod(sat);

plot_flag=[120 1 1 1 1 0 -4 3]; %   [samplerate considerNaN max(abs()) log uniform_colorbar power colorbar_limit1 colorbar_limit2]  
% plot_flag=[120 1 1 0 1 0 0 10]; %   [samplerate considerNaN max(abs()) log uniform_colorbar power colorbar_limit1 colorbar_limit2]  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEQCSPEC
% Purpose: short-session multipath processing and display
%          currently suitable for short-sessions (~2 hours) 
%          with data sampled at 1 sec <= rate <= 10 sec.
% 
% inputs:
%         time    = [t1 t2] - start and stop time in seconds from 00:00:00 UTC
%         max. plot freq. = [frqlim] - e.g. 0.2 Hz for ~0.05-0.06 Hz band [Default]
%         freq_cuttoff = [cutff1 cutff2] - freq. band of interest
%
% external functions: FONT FSA REPLACE TEQCPLOT POLARVIEW PLOTCLR BANDPASS 
%
% contact clement.ogaja@gmail.com 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectral analysis (Clement Ogaja, Geoscience Australia, 2006)
% References:
%            Ogaja & Satirapod, (2006). Analysis of high-frequency multipath in 1-Hz GPS 
%                                       kinematic solutions, submitted to GPS Solutions.
%            Ogaja & Hedfors, (2006). TEQC multipath metrics in MATLAB, submitted to GPS Toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file was last updated on February 22, 2009 to work with Matlab
% R2008b Release. The following functions were added: Checkfile.m, Scanfile.m,
% and Readfile.m

% History
% 03 Nov 2006 created using MATLAB 7.2 (R2006a)
% 16 Nov 2006 Enhanced to plot satellite tracks of multipath in a defined time range 
% 13 Dec 2006 Modified the code to remove system requirement of the
%             Signal Processing Toolbox
% 22 Feb 2009 Modified to work with R2008b. New functions added
%            (Checkfile.m, Scanfile.m, Readfiel.m)
% 06 Dec 2010 Modified to work with R2010a. Plotting control added
%              Glonass Support added


%
%
%%%%% Import multipath data  %%%%%%%%%%%%%%%%%%%%%%
%
if nargin==0
    [filen,path]=uigetfile('*.sn1;*.sn2;*.mp1;*.mp2',...
        'Pick your TEQC report file');
else
    [path,filen,ext]=fileparts(files);
    %path=[path '\'];
    fs=filesep;
    if ~isempty(path), path=[path fs]; end
    filen={[filen ext]};
end
file=char(filen);
%lkb
Constant.PlotName =[filen(1:end-4),'_']; %no extensions

global Sat_Capacity
Sat_Capacity=64;
global SatList
SatList=cell(1,Sat_Capacity);
for k=1:Sat_Capacity
    if k>=1 && k<=32
        SatList{k}=strcat('G',num2str(k,'%02d'));
    elseif k>=33 && k<=62
        SatList{k}=strcat('R',num2str(k-32,'%02d'));
    end
end
SatList{63}='E19'; %Galileo
SatList{64}='E20';

%%%%%%%%%%%%%%%%%%%%%%%%%%
N = checkfile([path file]);     % get totle line number
%%%%%%%%%%%%%%%%%%%%%%%%%%
if N<4 
    msgbox('invalid file -- too short!')
    return
end

n=0;i=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%
A = scanfile([path file],i);      % not really necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%

[t_samp mjl SAT n] = readfile(N,n,i,A,filen,[path file]);

T_SAMP=str2num(t_samp(max(find(t_samp==' ')):end));
MJL_START=str2num(mjl(max(find(mjl==' ')):end));
MJL(1)=MJL_START;
MJL(2:n)=MJL_START+(2:n)*T_SAMP/60/60/24;

JD=mjl2jd(MJL);
[type,maxy,miny]=get_filetype(file);

sat=SAT(1:n,:);
satval=sat; % multipath time series

%%%%% Import azimuth and elevation data  %%%%%%%%%%%%%%%%%%%%%%

out=teqcplot([path,(file(1:end-4)),'.azi']);az=out.azi;az_all=az;
	screen2print([Constant.PlotName 'Az'],Constant.PlotResolution); %save plots
out=teqcplot([path,(file(1:end-4)),'.ele']);el=out.ele;el_all=el;
	screen2print([Constant.PlotName 'El'],Constant.PlotResolution); %save plots
	
save midresult
load midresult




T1=1;
T2=min(172800,size(satval,1));
% T2=172796;
% cutff1=0.003;
% cutff2=0.016;
% cutff1=0.04;
% cutff2=0.05;
% plot_flag=[120 1 1 1 1 0 -4 3]; %   [samplerate considerNaN max(abs()) log uniform_colorbar power colorbar_limit1 colorbar_limit2]  
 

[t1 t2]=FindBestPeriod(sat);

%%%%% Extract data for satellites visible in the time window
vsats = 0;
% [tspan num_sats]=size(satval);
for i=1:Sat_Capacity,
    if (sum(isnan(satval(T1:T2,i)))/length(satval(T1:T2,i)))~=1,
		vsats=vsats+1;
		PRN(vsats)=i; % visible satellites
		MULT(:,vsats)=satval(T1:T2,i); % time window
        THETA(:,vsats)=az_all(T1:T2,i).*(pi/180); % time window
        RHO(:,vsats)=abs(el_all(T1:T2,i)-90)/90; % time window
        
        RHO(:,vsats)=RHO(:,vsats)-(RHO(:,vsats)>1).*(RHO(:,vsats)-1);
	end
end

statis

%%%%% Extract data for all satellites visible  
vsats_all = 0;
for i=1:Sat_Capacity,
    if (sum(isnan(satval(:,i)))/length(satval(:,i)))~=1,
		vsats_all=vsats_all+1;
		allPRN(vsats_all)=i; 
		allMULT(:,vsats_all)=satval(:,i); 
        for t=1:length(allMULT(:,vsats_all))
            JDn(t,vsats_all)=JD(t);
            if isnan(allMULT(t,vsats_all)), 
                JDn(t,vsats_all)=NaN; 
            end
        end
	end
end

%%%%% Compute FFT spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~frqlim,
    frqlim=(1/T_SAMP)*0.5;
end
for i=1:Sat_Capacity,
	Sout=fsa(satval(t1:t2,i),T_SAMP);
	freq(:,i)=Sout.f;
	Amp(:,i)=Sout.amp;
end

% Apply bandpass filter to extract high-frequency multipath in the band of
% interest
for i=1:vsats, %loop over all satellites
    MULT(:,i)=replace(MULT(:,i),NaN,0);
    HFMULT(:,i)=bandpass(MULT(:,i),cutff1,cutff2,length(MULT(:,i)),1);
    MULT(:,i)=replace(MULT(:,i),0,NaN);
    HFMULT(:,i)=replace(HFMULT(:,i),0,NaN);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% (1) display FFT spectra for all visible satellites,  %%%%%%%%%%%%%%%%%%%%%%
%%%%%     identify FFT peak(s) of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  %% raw time series ++++++++++++
map=colormap(flipud(jet));
miv=1;mav=Sat_Capacity;
%for j=1:vsats,
for j=vsats_all:-1:1,
    for i=1:Sat_Capacity,
        if i==allPRN(j),
            in=mod(round((i-1)*(length(map)-1)/(mav-miv)),8)*8+1;
            %--- Catch the out-of-range numbers
            if in==0;in=1;end
            if in > length(map);in=length(map);end
	        plot(JDn(:,j),allMULT(:,j)+allPRN(j),'color',map(in,:));hold on;
        end
    end
end
set(gca,'ylim',[0 Sat_Capacity+1])
set(gca,'ytick',[1:1:Sat_Capacity])
set(gca,'yticklabel',cell2mat(SatList'));
set(gca,'fontsize',7);
set(gca,'xticklabel',[JD]);
set(gca,'xlim',[JD(1) JD(end)])
datetick('x','HH:MM:SS','keeplimits','keepticks')
xlabel([datestr(JD(1)) '   |--------  T Samp: ' num2str(T_SAMP) ' s  --------|    ' datestr(JD(end))])
ylabel('SVs')
T=title(['TEQC Report file: ' strrep(file,'_','-')]);set(T,'fontsize',8)

screen2print([Constant.PlotName 'FFTraw'],Constant.PlotResolution); %save plots

figure;
h=axes('position',[0.11, 0.75, 0.85, 0.18]);  
map=colormap(flipud(jet));%colormap;
miv=1;mav=Sat_Capacity;
for i=1:Sat_Capacity,
    in=round((i-1)*(length(map)-1)/(mav-miv));
    %--- Catch the out-of-range numbers
    if in==0;in=1;end
    if in > length(map);in=length(map);end
	plot(freq(:,i),Amp(:,i),'color',map(in,:));hold on;
end
axis([0 frqlim 0 max(max(Amp))*1.25]);
% Re-format the colorbar
set(gca,'fontsize',7);
cb=colorbar('vertical');
%POS=get(cb,'position');
%set(cb,'position',[POS(1) POS(2) 0.03 POS(4)]);
set(cb,'fontsize',8);
set(get(cb,'xlabel'),'string','SVs');
set(cb,'box','on')
set(cb,'ylim',[1 length(map)]);
yal=linspace(1,length(map),4);
set(cb,'ytick',yal);
% Create the yticklabels
ytl=linspace(miv,mav,4);
s=char(4,4);
for i=1:4
    if min(abs(ytl)) >= 0.001
        B=sprintf('%2.0f',ytl(i));
    else
        B=sprintf('%2.0E',ytl(i));
    end
    s(i,1:length(B))=B;
end
set(cb,'yticklabel',s);
grid on
xlabel('Freq. (Hz)');
ylabel('Amplitude [m]');
T=title(['TEQC Report file: ' strrep(file,'_','-')]);set(T,'fontsize',8)
h=axes('position',[0.11, 0.11, 0.85, 0.55]); 
pcolor(freq',1:length(map),Amp');shading flat;cbar('v',[0 max(max(Amp))*1.25],'[m]');

[fn fu]=size(freq);
max_Amp1=max(max(Amp));max_Amp2=0;
for i=1:Sat_Capacity,
    freq(:,i)=replace(freq(:,i),NaN,0);
    Amp(:,i)=replace(Amp(:,i),NaN,0);
    
    if max(Amp(:,i)) >= max_Amp1,
	  f1=find(Amp(:,i)==max(Amp(:,i)));Frq1=freq(f1,i);max_sv1=i; max_Amp1=max(Amp(:,i));
    end


    for j=1:fn,
        if(freq(j,i)>=cutff1 && freq(j,i)<=cutff2),
            if(Amp(j,i)>=max_Amp2),max_Amp2=Amp(j,i);Frq2=freq(j,i);max_sv2=i;end;
        end
    end
end

 hold on;
 
 axis([0 frqlim 1 Sat_Capacity+1]);
 
% h=text(Frq1-0.005,max_sv1+.85,'o'); font(h,33);
% h=text(Frq2-0.003,max_sv2+.85,'o'); font(h,33); 
%to prevent crashing with really bad data
if exist('Frq1','var') 
   h=text(Frq1-0.002,max_sv1+.85,'o'); 
   font(h,33);  
end
if exist('Frq2','var') 
   h=text(Frq2-0.002,max_sv2+.85,'o'); 
   font(h,33);  
end
set(gca,'ylim',[1 Sat_Capacity+1])
set(gca,'ytick',[1.5:1:(Sat_Capacity+0.5)])
set(gca,'yticklabel',cell2mat(SatList'));
set(gca,'fontsize',7);

xlabel([datestr(JD(t1)) '   |--------  Freq. (Hz)  --------|    ' datestr(JD(t2))])
ylabel('SVS');
grid on

screen2print([Constant.PlotName 'FFTproc'],Constant.PlotResolution); %save plots


%%%%% (2) display unfiltered and filtered multipath  %%%%%%%%%%%%%%%%%%%%%%


temp=1:plot_flag(1):size(THETA,1);

temp_MULT=MULT(temp,:);
temp_HFMULT=HFMULT(temp,:);


if plot_flag(3)==1
    for k=1:length(temp)
%         temp_MULT(k,:)=max(abs(MULT(((k-1)*120+1):min([(k-1)*120+plot_flag(1) length(MULT)]),:)));
%         temp_HFMULT(k,:)=max(abs(HFMULT(((k-1)*120+1):min([(k-1)*120+plot_flag(1) length(HFMULT)]),:)));

        if plot_flag(1)==1
            temp_MULT(k,:)=abs(MULT(((k-1)*plot_flag(1)+1):min([(k-1)*plot_flag(1)+plot_flag(1) length(MULT)]),:));
            temp_HFMULT(k,:)=abs(HFMULT(((k-1)*plot_flag(1)+1):min([(k-1)*plot_flag(1)+plot_flag(1) length(HFMULT)]),:));
        else
            temp_MULT(k,:)=max(abs(MULT(((k-1)*plot_flag(1)+1):min([(k-1)*plot_flag(1)+plot_flag(1) length(MULT)]),:)));
            temp_HFMULT(k,:)=max(abs(HFMULT(((k-1)*plot_flag(1)+1):min([(k-1)*plot_flag(1)+plot_flag(1) length(HFMULT)]),:)));
        end
    end
    if plot_flag(4)==1
        temp_MULT=log10(temp_MULT);
        temp_HFMULT=log10(temp_HFMULT);
    end
    if plot_flag(6)==1
        temp_MULT=sqrt(temp_MULT);
        temp_HFMULT=sqrt(temp_HFMULT);
    end
end
if plot_flag(2)==1
    for k=1:length(temp)
        for kk=1:size(MULT,2)
            if isnan(temp_MULT(k,kk))
                temp_HFMULT(k,kk)=nan;
            end
        end
    end    
end

if plot_flag(5)==1
%     color_define=[-max(max(temp_MULT)) max(max(temp_MULT))];
    color_define=[plot_flag(7) plot_flag(8)];
end


figure
if plot_flag(5)==1
    polarview(THETA(temp,:),RHO(temp,:),temp_MULT,PRN,color_define);
    if plot_flag(4)==1
        cbar('v',color_define,'log(meter)')
    else
        cbar('v',color_define,'meter')
    end
else
    polarview(THETA(temp,:),RHO(temp,:),temp_MULT,PRN);
end
xlabel([datestr(JD(T1)) '   |--------  Unfiltered  --------|    ' datestr(JD(T2))])
T=title(['TEQC Report file: ' strrep(file,'_','-')]);set(T,'fontsize',8)

screen2print([Constant.PlotName 'MPraw'],Constant.PlotResolution); %save plot

figure;

if plot_flag(5)==1
    polarview(THETA(temp,:),RHO(temp,:),temp_HFMULT,PRN,color_define);

    if plot_flag(4)==1
        cbar('v',color_define,'log(meter)')
    else
        cbar('v',color_define,'meter')
    end
else
    polarview(THETA(temp,:),RHO(temp,:),temp_HFMULT,PRN);
end
xlabel([datestr(JD(T1)) '   |--------  Filtered  --------|    ' datestr(JD(T2))])
T=title(['TEQC Report file (debuged) : ' strrep(file,'_','-')]);set(T,'fontsize',8)

screen2print([Constant.PlotName 'MPproc'],Constant.PlotResolution); %save plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% FUNCTION: MJL DAYS TO JULIAN DAYS +++++++++++++++++++++++++++++++

function [out]=mjl2jd(in)

out=in+678941.999999741;

% FUNCTION: GET FILETYPE INFO +++++++++++++++++++++++++++++++++++++

function [out,maxy,miny]=get_filetype(teqfile);

% [path,name,ext,ver]=fileparts(teqfile);    %older version
[path,name,ext]=fileparts(teqfile);
switch ext
    case '.sn1'
        out='Signal to noise ratio S/N L1';
        maxy=10;
        miny=0;    
    case '.sn2'
        out='Signal to noise ratio S/N L2';
        maxy=10;
        miny=0;
    case '.mp1'
        out='Multipath L1';
        maxy=1;
        miny=-1;
    case '.mp2'
        out='Multipath L2';
        maxy=1;
        miny=-1;
    case '.iod'
        out='Derivative of ionospheric delay observable (m/s)';
        maxy=1;
        miny=-1;
    case '.ion'
        maxy=2;
        miny=-2;
        out='Ionospheric delay observable (m)';
    case '.ele'
        maxy=90;
        miny=0;
        out='Satellite elevation data';
    case '.azi'
        maxy=180;
        miny=-180;
        out='Satellite azimuthal data';
    otherwise
        disp('Somethings wrong..!')
end

% FUNCTION: PLACE A MODIFIED COLORBAR ++++++++++++++++++++++++++++++

function CB=cbar(loc,range,label);

% .............................................................
% CB = cbar(loc,range,label)
%   places a colorbar at:
%   loc = 'v' in vertical or 'h' in horizontal
%           position in current figure scaled between:
%   range = [min max] with a:
%   label = 'string'.
%
%   fontsize is reduced to 10 and width of bar is half default.
%   
%   Example:    [X,Y,Z]=peaks(25);
%               range=[min(min(Z)) max(max(Z))];
%               pcolor(X,Y,Z);
%               cbar('v',range,'Elevation (m)')
% .............................................................

caxis([range(1) range(2)]);
switch loc
    case 'v'
        CB=colorbar('vertical');
        set(CB,'ylim',[range(1) range(2)]);
        %POS=get(CB,'position');
        %set(CB,'position',[POS(1) POS(2) 0.03 POS(4)]);
        set(CB,'fontsize',8);
        set(get(CB,'xlabel'),'string',label);
        set(CB,'box','on')

    case 'h'
        CB=colorbar('horizontal');
        set(CB,'xlim',[range(1) range(2)]);
        %POS=get(CB,'position');
        %set(CB,'position',[POS(1) POS(2) POS(3) 0.03]);
        set(CB,'fontsize',8);
        set(get(CB,'xlabel'),'string',label)
end

% FUNCTION: SECONDS TO HOURS, MINUTES and SECONDS ++++++++++++++++

function timestr=secs2hms(SECS)

HOURS=SECS/60/60;
hours=floor(HOURS);
MINUTES=(HOURS-hours)*60;
minutes=floor(MINUTES);
seconds=(MINUTES-minutes)*60;
HH=num2str(hours);
MM=num2str(minutes);
SS=num2str(seconds);

if seconds<10;
    SS=['0' num2str(SS)];
else
    SS=num2str(SS);
end
if minutes<10;
    MM=['0' num2str(MM)];
else
    MM=num2str(MM);
end
timestr=[HH ':' MM ':' SS];

% EOF +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



