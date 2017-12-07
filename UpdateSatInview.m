function A=UpdateSatInview(tline)
    % note: for readfile only, not for readfile_v2

    global SatList
    
    temp=strread(tline,'%s');
    
    
    if isempty(temp) 
        A=[];
    else
        x=str2num(temp{1});
        switch x
            case {0,-1}
                A=x;
            otherwise
                if length(temp)==x+1
                    A=zeros(1,x);
                    for k=2:(x+1)
                        A(k-1)=find(strcmp(SatList,temp(k)));
                    end
                else
                    A=[];
                end
        end
        
    end
        
end