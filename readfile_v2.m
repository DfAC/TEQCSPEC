function [t_samp,TimeStamp,SatVal] = readfile_v2(N,num_Headerlines,file)
% to read compact3 data file from TEQC 
%
% See also TEQCSPEC_main, checkfile
%
% History
% 22 Feb 2009 created using Matlab R2008b
% 07 Nov 2017 major modified using Matlab R2014b by Lei Yang 

    global Sat_Capacity     % 1-32:GPS; 33-64:GLONASS; 65-96:GALILEO; 97-99:SBAS; 100:Reserved; 101-140:BeiDou; 141-150:QZSS; 151-160:IRNSS

    SatVal=nan(N,Sat_Capacity);  
    TimeStamp=zeros(1,N);
    Ind_line=0;                    % index of line
    Ind_epoch=0;                    % index of epoch
    

    fid=fopen(file);
    % h=waitbar(0,['Please wait... Importing ' file]);
    
    while ~feof(fid);
        tline = fgetl(fid);     % get header or SV info line
        Ind_line=Ind_line+1;
        
        if Ind_line<=num_Headerlines
            temp=textscan(tline,'%s',1);
            if strcmp(temp{1},'GPS_START_TIME')     % take the information for the starting time
                 Line_Data=textscan(tline,'%s %n %n %n %n %n %n',1);
                 TimeStamp=TimeStamp+datenum([Line_Data{2:7}])*86400;
            end           
        else
    %         waitbar(Ind_line/N,h);drawnow    
            if mod(Ind_line,5000)==0
                disp(['Reading Data -- ' strcat(num2str(100*Ind_line/N),'%')])
            end
            
            Ind_epoch=Ind_epoch+1;
            
            % update current timestamp
            [temp,pos]=textscan(tline,'%n',1);
            if Ind_epoch==1
            elseif Ind_epoch==2
                t_samp=temp{1}-time_step;           % take the first interval as the official interval
                t_samp=round(t_samp*100)/100;       % round up to 10 milli-sec
            else
                if abs((temp{1}-time_step)-t_samp)>0.01
%                     dbstop at 47
                    disp(['Warning: Time stamp interval mismatch! check ' file ' line ' num2str(Ind_line)]);
                end
            end
            time_step=temp{1};
            TimeStamp(Ind_epoch)=TimeStamp(Ind_epoch)+time_step;
            
            
            if cell2mat(textscan(tline(pos+1:end),'%n',1))~=-1
                ActiveSV=UpdateSatInview(tline(pos+1:end));
            end
            
            tline = fgetl(fid); % get the data line
            Ind_line=Ind_line+1;
            
            SatVal(Ind_epoch,ActiveSV)=cell2mat(textscan(tline,'%n',length(ActiveSV)))';    % record data
            
        end
    end
    
    fclose(fid);
    % close(h);
    
    SatVal=SatVal(1:Ind_epoch,:);
    TimeStamp=TimeStamp(1:Ind_epoch);
end

function A=UpdateSatInview(tline)

    global SatList
    
    [temp, pos]=textscan(tline,'%n',1);
    if isempty(temp), A=[]; return;    end
    
    ID_SV=textscan(tline(pos+1:end),'%s',temp{1});
    
    A=cellfun(@(s) find(strcmp(SatList, s)), ID_SV{1})';
        
end

