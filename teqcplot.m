function out=teqcplot(files,T_SAMP,TimeStamp,SatVal,type,maxy,miny)

% TEQCPLOT  Read and plot TEQC report files
%
%   TEQC is the Toolkit for GPS/GLONASS/GALILEO/BeiDou/SBAS/QZSS/IRNSS Data used to solve 
%   many pre-processing problems with GPS, GLONASS, and SBAS data:
%   TEQC stands for Translation, Editing, and Quality Check
%   More info here: http://facility.unavco.org/software/teqc/teqc.html
%
%   TEQCPLOT(FILENAME) returns the content of a TEQC reportfile in
%   a figure (and a struct variable).
%   Input argument FILENAME is optional. 
%
%   Valid TEQC report files are:
%  
%   FORMAT  DESCRIPTION
%   *.sn1   Signal to noise ratio (S/N) Carrier L1
%   *.sn2   Signal to noise ratio (S/N) Carrier L2
%   *.iod   Derivative of ionospheric delay observable (m/s)
%   *.ion   Ionospheric delay observable (m)
%   *.mp1   Multipath Carrier L1
%   *.mp2   Multipath Carrier L2
%   *.azi   Satellite azimuthal data (degrees)
%   *.ele   Satellite elevation data (degrees)
%
%   Jim Hedfors, Uppsala University, 2005
%   Comments and requests for more TEQC Matlab tools to:
%   jim.hedfors@gmail.com

% FUNCTION: READER & PLOTTER ++++++++++++++++++++++++++++++++++++++++++++++

    if nargin==0
        [temp,path]=uigetfile('*.sn1;*.sn2;*.iod;*.ion;*.mp1;*.mp2;*.m12;*.m15;*.m21;*.mp1;*.mp2;*.azi;*.ele',...
            'Pick your TEQC report file');

        [~, file, ext]=fileparts(temp);
        clear temp
    else
        [path,file,ext]=fileparts(files);
        if isempty(ext), disp('Please include file extension.'); return ; end
        if isempty(path), path=pwd; end            
    end
    if path(end)~=filesep
        path=[path filesep];
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin<=1

        [N, file_version] = checkfile([path file ext]);     % get totle line number
        if N<4 
            msgbox('invalid file -- too short!')
            return
        end

        switch file_version{1}
            case 'COMPACT3'
                num_Headerlines=2;
                [T_SAMP TimeStamp SatVal] = readfile_v2(N,num_Headerlines,[path file ext]);
                % T_SAMP -- epoch time interval; TimeStamp -- epoch timestamp list; SatVal -- main data for each SV ; 
            otherwise       % for older datafile version
                n=0;i=4;
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                A = scanfile([path file],i);      % not really necessary
                %%%%%%%%%%%%%%%%%%%%%%%%%%

                [t_samp mjl SatVal n] = readfile(N,n,i,A,[path file]);

                T_SAMP=str2num(t_samp(max(find(t_samp==' ')):end));
                TimeStamp=mjl2date(str2num(mjl(max(find(mjl==' ')):end)))*86400+(1:n)*T_SAMP;
                SatVal=SatVal(1:n,:);      % multipath time series
        end
    end
    if nargin<=4
        [type,maxy,miny]=get_filetype([path file ext]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    ActiveSV_ind=find(sum(~isnan(SatVal)));

    global Sat_Capacity SatList
    figure;box on;hold on
    pcolor(TimeStamp',1:length(ActiveSV_ind),SatVal(:,ActiveSV_ind)')   % note: the unit of x axis is actually in second
    set(gcf,'renderer','zbuffer');
    shading flat
    
    if TimeStamp(end)-TimeStamp(1)<6*3600
        xtick_interval=3600;   % unit : sec
    elseif TimeStamp(end)-TimeStamp(1)<24*3600
        xtick_interval=3*3600;   % unit : sec
    else
        xtick_interval=12*3600;  % unit : sec
    end
    temp=(ceil(TimeStamp(1)/xtick_interval):floor(TimeStamp(end)/xtick_interval))*xtick_interval;
    set(gca,'xtick',temp)
    set(gca,'xticklabel',datestr(temp/86400,'HH:MM'));
    clear temp
    set(gca,'xlim',[TimeStamp(1) TimeStamp(end)])
    set(gca,'ylim',[1 length(ActiveSV_ind)+1])
    set(gca,'ytick',[1.5:1:(Sat_Capacity+0.5)])
    set(gca,'yticklabel',cell2mat(SatList(ActiveSV_ind)'));
    set(gca,'fontsize',7);
    xlabel([datestr(TimeStamp(1)/86400) '   |--------  T Samp: ' num2str(T_SAMP) ' s  --------|    ' datestr(TimeStamp(end)/86400)])
    ylabel('SVs')
    T=title(['TEQC Report file: ' strrep(file,'_','-')]);
    set(T,'fontsize',8)

    cbar('v',[miny maxy],type);
    colormap(flipud(jet));
    caxis([miny maxy]);

    out.(ext(end-2:end))=SatVal;
    out.T_samp=T_SAMP;
    out.TimeStamp=TimeStamp;
end

% FUNCTION: MJL DAYS TO Date  +++++++++++++++++++++++++++++++
function [out]=mjl2date(in)
    out=in + datenum([1858 11 17]);
end

% FUNCTION: GET FILETYPE INFO +++++++++++++++++++++++++++++++++++++
function [out,maxy,miny]=get_filetype(teqfile)

    % [path,name,ext,ver]=fileparts(teqfile);       % old version
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
            out='Satellite Elevation';
        case '.azi'
            maxy=180;
            miny=-180;
            out='Satellite Azimuthal';
        otherwise
            disp('Somethings wrong..!')
    end
end

% FUNCTION: PLACE A MODIFIED COLORBAR ++++++++++++++++++++++++++++++
function CB=cbar(loc,range,label)

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
            set(get(CB,'ylabel'),'string',label);
            set(CB,'box','on')

        case 'h'
            CB=colorbar('horizontal');
            set(CB,'xlim',[range(1) range(2)]);
            %POS=get(CB,'position');
            %set(CB,'position',[POS(1) POS(2) POS(3) 0.03]);
            set(CB,'fontsize',8);
            set(get(CB,'xlabel'),'string',label)
    end
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
end

% EOF +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++