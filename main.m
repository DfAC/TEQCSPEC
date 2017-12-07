function out=TEQCSPEC_main(files);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   TEQCSPEC  TEQC multipath spectrum analysis
    %   ++++++++++++++++++++++++++++++++
    %   TEQCSPEC(FILENAME) makes a geographically oriented plot of 
    %   multipath/SNR in a re-configureable band 0.005 - 0.008 Hz (default) in a defined time range 
    %
    %   History
    %   03 Nov 2006 created using MATLAB 7.2 (R2006a) by clement.ogaja@gmail.com
    %   16 Nov 2006 Enhanced to plot satellite tracks of multipath in a defined time range 
    %   13 Dec 2006 Modified the code to remove system requirement of the Signal Processing Toolbox
    %
    %   22 Feb 2009 Modified to work with R2008b. New functions added (Checkfile.m, Scanfile.m, Readfiel.m)
    %   06 Dec 2010 Modified to work with R2010a. Plotting control added, Glonass Support added
    %   2010 -2015 maintained by Lei Yang and Bonenberg Lukasz 
    %
    %   2017 an major update by Lei Yang
    %     -- support latest TEQC (version 15/09/2017) and Matlab R2014b
    %     -- support COMPACT3 input data (multipaht and SN)
    %     -- support GPS/GLONASS/GALILEO/BeiDou/SBAS/QZSS/IRNSS
    %     -- the bandpass filter could be applied to just the multipath SVs
    %     -- improved plotting function
    %     -- centralized control area
    %     -- more comments added
    %     -- scripts house keeping
    %
    %   plot for a period up to 48hr
    %
    %   ++++++++++++++++++++++++++++++++
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  +++++++++++++++++++++++++++++++++++
    %%%  control segment
    %%%%%%%%%%%%%%%%%%%%%%
    
    close(findall(0,'Type','figure')); 
%     clear all
      
    % 
    dbstop if error                   % for debugging
%     dbstop at 146                     % for break after spectrum analysis, remember to change this after each update
    
    run_data_reading=1;                 % use this indicator if the midresult is already generated from a earlier run, pure to save time
        workspace_afterdatareading='midresult_afterdatareading.mat';    % will load this file if run_data_reading==0; similar for below
%         workspace_afterdatareading='midresult_1913.mat';     
    run_spectrum_analysis=1 ;   
        workspace_afterspecanlz='midresult_afterspecanlz';              % will load this file if run_spectrum_analysis==0; similar for below
    plot_skyplot=1;
    save_plot=1;

    %%%%***********************
    % set up start/end time window, or use full period (up to first 48hr) by setting nan value
    %
    % NOTE(*IMPORTANT*): 
    %
    % To identify the true multipath, you need to
    % 1) run full period investigation, to find the suspectious multipath time slot
    % 2) change the below start/end pair (unit:epoch), and run the scripts for 2nd time, to investigate
    % the multipath period only, and only at this time the frequency analysis is meaningful 
    % 3) apply the Bandpass filter and run the scripts for 3rd time to surpress the multipath freq, 
    % all sv at all time will be impacted, however the multipath period shall be impacted more effectively 
    % 4) alternatively, the bandpass filter could be applied only to the SVs defined in the SVind_BP, 
    % the index here is refered to SatList index
    T1=nan;                             
    T2=nan;                  % minimum length 200 epochs (as the freq sampling interval =1/Ts/(T2-T1+1) ), suggest 300 epochs ; the freq sampling interval used in plotting is 0.005Hz          
    SVind_BP=[];
    %%%%***********************
    
    
    
    % take the input file name
    if nargin==0
        [temp,path]=uigetfile('*.sn1;*.sn2;*.mp1;*.mp2;*.m12;*.m15;*.m21',...
            'Pick your TEQC report file for skyplot preview');
        %note that only selected file get plotted so you might be forced to run software 4 times

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
    PlotName_prefix =[file,'_']; %no extensions
    
    
    % freq setting
    frqlim = 0.2;                       % Max plotting frequency (Hz)
    frqstep=0.005;
    frq_level=frqlim/frqstep;             
    % logic : max freq is 1/2/Ts, equal to 0.5 if Ts=1Hz; the minimal
    % investigate period is 100 epochs (), thus the 
    cutoff_L = 0.005;                     % Lower cut off frequency (Hz)
    cutoff_H = 0.008;                     % Upper cut off frequency (Hz)%
    
    shortestperiod_freq_analyze=200;    % unit: epoch, the period shorter than this will not be put into freq analyzer

    
    % Plot setting
    Plot_Resolution=300;
    Plot_SampleRate=60;     % if the data size is large, plotting will be quite slow, and it is not necessary to plot all points,
                            % so 1 points will be plotted from every Plot_Plot_SampleRate epoches 
    Plot_SampleMethod=1;    % 0: the sample point value; 
                            % 1: max(abs(*)) of the sampling period; 
                            % 2: mean(abs(*)) of the sampling period
                            % 3: mean(*) of the sampling period
    Plot_LogValue=0;        % 0: normal value; 1: log10 value
    Plot_SqrtValue=0;       % 0: normal value; 1: square root value
    Plot_ValueRange=[0 10]; % out of range values will be colored as the value of range ends
    if ~isempty(strfind(ext,'.sn'))     % SNR default settings
        Plot_SampleMethod=2;
        Plot_ValueRange=[12 54];        % adjustable, default [12 54] (dBHz) according to the RTCM
    elseif ~isempty(strfind(ext,'.m'))  % Multipath default settings
        Plot_ValueRange=[1 0];
    end
        
    %%%%%%%%%%%%%%%%%%%%%%
    %%%  control segment end
    %%%  +++++++++++++++++++++++++++++++++++
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %% set up constant
    global Sat_Capacity
    Sat_Capacity=160;
    global SatList
    SatList=cell(1,Sat_Capacity);       % SV ID 
    for k=1:Sat_Capacity
        if k>=1 && k<=32
            SatList{k}=strcat('G',num2str(k,'%02d'));
        elseif k>=33 && k<=64
            SatList{k}=strcat('R',num2str(k-32,'%02d'));
        elseif k>=65 && k<=96
            SatList{k}=strcat('E',num2str(k-64,'%02d'));
        elseif k>=97 && k<=99
            switch k
                case 97
                    SatList{k}='S20';
                case 98
                    SatList{k}='S23';
                case 99
                    SatList{k}='S36';
            end
        elseif k==100
            % reserved 
            SatList{k}='rsv';
        elseif  k>=101 && k<=140
            SatList{k}=strcat('C',num2str(k-100,'%02d'));
        elseif  k>=141 && k<=150
            SatList{k}=strcat('J',num2str(k-140,'%02d'));
        elseif  k>=151 && k<=160
            SatList{k}=strcat('I',num2str(k-150,'%02d'));
        end
    end

    %%
    %%%%% Import multipath/sn data  %%%%%%%%%%%%%%%%%%%%%%
    if run_data_reading

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        [N, file_version] = checkfile([path file ext]);     % get total line number and version
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        if N<4 
            msgbox('invalid file -- too short!')
            return
        end

        switch file_version{1}
            case 'COMPACT3'
                num_Headerlines=2;
                [T_SAMP,TimeStamp,SatVal] = readfile_v2(N,num_Headerlines,[path file ext]);
                % T_SAMP -- epoch time interval; TimeStamp -- epoch timestamp list; SatVal -- main data for each SV ; 

            otherwise       % to support older datafile version (<COMPACT3), remove later
                n=0;i=4;
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                A = scanfile([path file],i);      % not really necessary
                %%%%%%%%%%%%%%%%%%%%%%%%%%

                [t_samp,mjl,SatVal,n] = readfile(N,n,i,A,[path file]);

                T_SAMP=str2num(t_samp(max(find(t_samp==' ')):end));
                TimeStamp=mjl2date(str2num(mjl(max(find(mjl==' ')):end)))*86400+(1:n)*T_SAMP;
                SatVal=SatVal(1:n,:);      % multipath time series
        end


    %     [type,maxy,miny]=get_filetype(file);      % get data range depending on the file type


        %%%%% Import azimuth and elevation data  %%%%%%%%%%%%%%%%%%%%%%

        out=teqcplot([path file '.azi']);az_all=out.azi;clear out; az_all=mod(az_all,360);
            if save_plot, screen2print([PlotName_prefix 'Azi'],Plot_Resolution); end   %save plots
        out=teqcplot([path file '.ele']);el_all=out.ele;clear out; 
            if save_plot, screen2print([PlotName_prefix 'Ele'],Plot_Resolution); end   %save plots
            
        save(workspace_afterdatareading,'file_version','T_SAMP','TimeStamp','SatVal','az_all','el_all')
    else
        load(workspace_afterdatareading)
    end

    %%
    if run_spectrum_analysis 

        if isnan(T1)    T1=1; end
        if isnan(T2)    T2=min(172800,size(SatVal,1)); end

%         [t1, t2]=FindBestPeriod(SatVal,SVinview_least); % multipath period to be analyzed, require at least [SVinview_least] sv in view

        %%%%% Extract data for satellites visible in the time window
        SV_visiable=find(sum(isnan(SatVal(T1:T2,:)))<size(SatVal(T1:T2,:),1));
        vsats = length(SV_visiable);
        % apply time window
        SatVal_Active=SatVal(T1:T2,SV_visiable); 
        THETA=az_all(T1:T2,SV_visiable).*(pi/180); 
        RHO=abs(el_all(T1:T2,SV_visiable)-90)/90; 
        RHO=RHO-(RHO>1).*(RHO-1);

        statis              % output basic statistics

        %%%%% Extract data for all satellites visible  (???)
        SV_visiable_fullperiod=find(sum(isnan(SatVal))<size(SatVal,1));
        vsats_fullperiod = length(SV_visiable_fullperiod);
        SatVal_Active_fullperiod=SatVal(:,SV_visiable_fullperiod);
        TimeStamp_withNaN=TimeStamp'*ones(1,vsats_fullperiod);
        TimeStamp_withNaN(find(isnan(SatVal_Active_fullperiod)))=nan;
        

        %%%%% Compute FFT spectra
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~frqlim,     frqlim=(1/T_SAMP)*0.5;    end
        freq=cell(1,Sat_Capacity);
        Amp=cell(1,Sat_Capacity);
        for i=SV_visiable
            [tempt1, tempt2]=FindBestPeriod(SatVal(:,i));

            if tempt2-tempt1>shortestperiod_freq_analyze
                Sout=fsa(SatVal(tempt1:tempt2,i),T_SAMP);
                freq{i}=Sout.f;
                Amp{i}=Sout.amp;
            end
        end
%         for i=1:Sat_Capacity
%             if sum(~isnan(SatVal(t1:t2,i)))>0
%                 [tempt1, tempt2]=FindBestPeriod(SatVal(t1:t2,i));
%                 
%                 if tempt2-tempt1>shortestperiod_freq_analyze
%                     Sout=fsa(SatVal(tempt1:tempt2,i),T_SAMP);
%                     freq{i}=Sout.f;
%                     Amp{i}=Sout.amp;
%                 end
%             end
%         end

        % Apply bandpass filter to extract high-frequency multipath in the band of interest
        if ~isempty(strfind(ext,'.m'))
            SatVal_Active_BP=nan(size(SatVal_Active));
            for i=1:vsats, %loop over all satellites
                temp_ind=find(~isnan(SatVal_Active(:,i)));

                j=1;
                while j+shortestperiod_freq_analyze<=length(temp_ind)
                    j_end=j+shortestperiod_freq_analyze;
                    while j_end<length(temp_ind) && temp_ind(j_end)-temp_ind(j)==j_end-j
                        j_end=j_end+1;
                    end
                    if j_end-1-j>shortestperiod_freq_analyze
                        % Note here: check the bandpass function, may be wrong!
                        % <<<<<<<<<<<<<<---------------------!!!!!!!!
                        SatVal_Active_BP(temp_ind(j:(j_end-1)),i)=bandpass(SatVal_Active(temp_ind(j:(j_end-1)),i),cutoff_L,cutoff_H,length(SatVal_Active(temp_ind(j:(j_end-1)),i)),1);
                        j=j_end;
                    else
                        j=j+1;
                    end
                end
            end
        end
%             SatVal_Active(:,i)=replace(SatVal_Active(:,i),NaN,0);
%             SatVal_Active_BP(:,i)=bandpass(SatVal_Active(:,i),cutoff_L,cutoff_H,length(SatVal_Active(:,i)),1);
%             SatVal_Active(:,i)=replace(SatVal_Active(:,i),0,NaN);
%             SatVal_Active_BP(:,i)=replace(SatVal_Active_BP(:,i),0,NaN);
    
    
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%  PLOTS
        %%%%% (1) display FFT spectra for all visible satellites,  %%%%%%%%%%%%%%%%%%%%%%
        %%%%%     identify FFT peak(s) of interest
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% raw time series ++++++++++++
        figure                 
        map_level=32;
        map=colormap(jet(map_level));
        
        for j=fliplr(SV_visiable_fullperiod)
            if j<=32
                subplot(2,2,1)
            elseif j<=64
                subplot(2,2,2)
            elseif j<=96
                subplot(2,2,3)
            else
                subplot(2,2,4)
            end
            
            hold on
            
            in=mod(mod(j,map_level),16)*map_level/16+1;    
            i=find(SV_visiable_fullperiod==j);
            plot(TimeStamp_withNaN(:,i), ...
                (SatVal_Active_fullperiod(:,i)-mean(Plot_ValueRange))/(Plot_ValueRange(2)-Plot_ValueRange(1)) ...
                + mod(SV_visiable_fullperiod(i)-1,32)+1+ 128*(SV_visiable_fullperiod(i)>128),'color',map(in,:));
        end
        
        for j=1:4
            subplot(2,2,j)
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
            set(gca,'fontsize',7);
            xlabel([datestr(TimeStamp(1)/86400) '   |--------  T Samp: ' num2str(T_SAMP) ' s  --------|    ' datestr(TimeStamp(end)/86400)])
            
            switch j
                case 1
                    set(gca,'yticklabel',cell2mat(SatList(1:32)'));
                    set(gca,'ylim',[0 32+1])
                    set(gca,'ytick',[0.5:1:32])
                    
                    title('GPS','fontsize',8)
                case 2
                    set(gca,'yticklabel',cell2mat(SatList(33:64)'));
                    set(gca,'ylim',[0 32+1])
                    set(gca,'ytick',[0.5:1:32])
                    
                    title('GLONASS','fontsize',8)
                case 3
                    set(gca,'yticklabel',cell2mat(SatList(65:96)'));
                    set(gca,'ylim',[0 32+1])
                    set(gca,'ytick',[0.5:1:32])
                    
                    title('GALILEO','fontsize',8)
                case 4
                    set(gca,'yticklabel',cell2mat(SatList(97:160)'));
                    set(gca,'ylim',[0 64+1])
                    set(gca,'ytick',[0.5:1:64])
                    
                    title('BEIDOU/GEO/QZSS/IRNSS','fontsize',8)
            end
        end
        subplot(2,2,1)
        text(TimeStamp(end),36,['TEQC Report file: ' strrep(file,'_','-')],'fontsize',12,'FontWeight','bold')
        
        if save_plot
            screen2print([PlotName_prefix 'FFTraw'],Plot_Resolution); %save plots
        end
                
        %% frequence plots +++++++++++++ multipath only
        if ~isempty(strfind(ext,'.m')) 
            figure;     
            h=axes('position',[0.1, 0.75, 0.85, 0.2]);  hold on;
            map_level=32;
            map=colormap(jet(map_level));
            for j=SV_visiable,
                in=mod(mod(j,map_level),16)*map_level/16+1;   
                %--- Catch the out-of-range numbers
                if in==0;in=1;end
                if in > length(map);in=length(map);end
                plot(freq{j},Amp{j},'color',map(in,:)); 
            end
            axis([0 frqlim 0 max(cell2mat(cellfun(@(x) max(x),Amp,'UniformOutput',false)))]);
            grid on
            xlabel('Freq. (Hz)');
            ylabel('Amplitude [m]');
            set(gca,'fontsize',7);
            T=title(['TEQC Report file: ' strrep(file,'_','-')]);set(T,'fontsize',8)

            h=axes('position',[0.1, 0.1, 0.85, 0.55]); 

            temp_freq=frqstep:frqstep:frqlim;
            temp_amp=nan(length(SV_visiable),frq_level);
            for j=SV_visiable
                if ~isempty(freq{j})
                    temp_amp(find(SV_visiable==j),:)=interp1(freq{j},Amp{j},temp_freq);
                end
            end
            temp_ind=find(sum(~isnan(temp_amp'))>0);


            surf(temp_freq-frqstep/2,(1:(length(temp_ind)+1))-0.5,temp_amp([temp_ind temp_ind(end)],:)); shading flat; view(2) % add the last raw to suit the need of 'surf'
            caxis([0 1])


            set(h,'ytick',1:(length(temp_ind)))
            set(h,'yticklabel',cell2mat(SatList(SV_visiable(temp_ind))'))
            set(h,'FontSize',7)
            xlabel('Freq. (Hz)');
            ylabel('Visible SVs');
            axis([0 frqlim 0 length(temp_ind)+1]);


            %% find the suspecious multipath sv and period 
            temp_Amp=cell(size(Amp));
            for i=1:Sat_Capacity,
                if isempty(freq{i}), continue; end

                t1=sum(freq{i}<=0.1*frqlim)+1; t2=sum(freq{i}<=frqlim); % exclude the 10% lower region to avoid noise
                temp_Amp{i}=Amp{i}(t1:t2);
            end
            [max_Amp max_Amp_row]=max(cell2mat(cellfun(@(x) max(x),temp_Amp,'UniformOutput',false)));
            temp=find(cell2mat(cellfun(@(x) ~isempty(x),Amp,'UniformOutput',false)));
            max_Amp_row=temp(max_Amp_row); clear temp
            max_Amp_column=find(temp_Amp{max_Amp_row},max_Amp);

            disp(['Max value in the range [ ', num2str(0.1*frqlim), ' ~ ' , num2str(frqlim),' ] is ', num2sr(max_Amp), ' on ', SatList{max_Amp_row}, ' @ ' num2str(freq{max_Amp_row}(max_Amp_column)) ' Hz'])


            %% save figure
            if save_plot
                screen2print([PlotName_prefix 'FFTproc'],Plot_Resolution); %save plots
            end
        end
        
        save(workspace_afterspecanlz)
    else
        load(workspace_afterspecanlz)
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  PLOTS
    %%%%% (2) display unfiltered and filtered multipath 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    sample_SatVal_Active=apply_SampleMode(SatVal_Active, Plot_SampleRate, Plot_SampleMethod,Plot_LogValue,Plot_SqrtValue);
    if ~isempty(strfind(ext,'.m'))
        sample_SatVal_Active_BP=apply_SampleMode(SatVal_Active_BP, Plot_SampleRate, Plot_SampleMethod,Plot_LogValue,Plot_SqrtValue);
        sample_SatVal_Active_BP(isnan(sample_SatVal_Active))=nan;
    end
    if Plot_SampleRate>1
        THETA=THETA(1:Plot_SampleRate:end,:);
        RHO=RHO(1:Plot_SampleRate:end,:);
    end
      
    if plot_skyplot
        plot_sky(sample_SatVal_Active,THETA,RHO,SV_visiable,Plot_ValueRange,Plot_LogValue,TimeStamp([T1,T2]),file,ext);
        xlabel([datestr(TimeStamp(T1)/86400) '   |--------  Unfiltered  --------|    ' datestr(TimeStamp(T2)/86400)])

        if save_plot
            screen2print([PlotName_prefix 'MPraw'],Plot_Resolution); %save plot
        end

        if ~isempty(strfind(ext,'.m')) % only multipath plot has a filtered version
            plot_sky(sample_SatVal_Active_BP,THETA,RHO,SV_visiable,Plot_ValueRange,Plot_LogValue,TimeStamp([T1,T2]),file,ext);
            xlabel([datestr(TimeStamp(T1)/86400) '   |--------  Filtered  --------|    ' datestr(TimeStamp(T2)/86400)])
            if save_plot
                screen2print([PlotName_prefix 'MPproc'],Plot_Resolution); %save plot
            end
        end
    end
end


% FUNCTION: MJL DAYS TO Date  +++++++++++++++++++++++++++++++
function [out]=mjl2date(in)
    out=in + datenum([1858 11 17]);
end

% FUNCTION: GET FILETYPE INFO +++++++++++++++++++++++++++++++++++++
function [out,maxy,miny]=get_filetype(teqfile);

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
end

% FUNCTION: PLACE A MODIFIED COLORBAR +++++++++++++++++++++++++++++
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

    caxis([min(range) max(range)]);
    switch loc
        case 'v'
            CB=colorbar('vertical');
            set(CB,'ylim',[min(range) max(range)]);
            %POS=get(CB,'position');
            %set(CB,'position',[POS(1) POS(2) 0.03 POS(4)]);
            set(CB,'fontsize',8);
            set(get(CB,'xlabel'),'string',label);
            set(CB,'box','on')

        case 'h'
            CB=colorbar('horizontal');
            set(CB,'xlim',[min(range) max(range)]);
            %POS=get(CB,'position');
            %set(CB,'position',[POS(1) POS(2) POS(3) 0.03]);
            set(CB,'fontsize',8);
            set(get(CB,'xlabel'),'string',label)
    end
end

% FUNCTION: SECONDS TO HOURS, MINUTES and SECONDS +++++++++++++++
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

function output=apply_SampleMode(input,Plot_SampleRate, Plot_SampleMethod,Plot_LogValue,Plot_SqrtValue)

    if Plot_SampleMethod==1 || Plot_SampleMethod==2,    input=abs(input);    end

    if Plot_SampleRate>1 
        if Plot_SampleMethod>=1 
            while mod(size(input,1),Plot_SampleRate)>0
                input=[input; input(end,:)];
            end
            
            temp=reshape(input,Plot_SampleRate, size(input,1)/Plot_SampleRate, size(input,2));
            if Plot_SampleMethod==1
                temp=max(temp);
            else
                temp=mean(temp);
            end
            output=reshape(temp,size(input,1)/Plot_SampleRate,size(input,2));
            
%             for k=1:ceil(size(input,1)/Plot_SampleRate)
%                 input_slide=input(((k-1)*Plot_SampleRate+1):min([k*Plot_SampleRate length(input)]),:);
% 
%                 if Plot_SampleMethod==1
%                     output(k,:)=max(input_slide);
%                 else    % Plot_SampleMethod==2 || Plot_SampleMethod==3
%                     output(k,:)=mean(input_slide);
%                 end
%             end
%             clear input_slide 
        else
            output=input(1:Plot_SampleRate:end,:);
        end
    else
        output=input;
    end
    
    if Plot_LogValue==1,    output=log10(output);   end
    if Plot_SqrtValue==1,   output=sqrt(output);    end
end

function plot_sky(sample_SatVal_Active,THETA,RHO,SV_visiable,Plot_ValueRange,Plot_LogValue,TimeStamp,file,ext)
        
        figure
        
        polarview(THETA,RHO,sample_SatVal_Active,SV_visiable,Plot_ValueRange);
        
        set(get(gca,'ylabel'),'string','Azimuth (deg)');

 
        switch ext(1:3)
            case '.sn'
                plot_unit='unitless';
            case '.mp'
                plot_unit='meter';
            case '.io'
                plot_unit='meter';
            otherwise
                disp('Somethings wrong!')
        end

        if Plot_LogValue==1
            cbar('v',Plot_ValueRange,['log(' plot_unit ')'])
        else
            cbar('v',Plot_ValueRange,plot_unit)
        end
        
        T=title(['TEQC Report file: ' strrep(file,'_','-')]);set(T,'fontsize',8)
end

% EOF ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



