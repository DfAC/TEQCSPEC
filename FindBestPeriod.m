function [t1,t2]=FindBestPeriod(sats,SVinview_least)

if nargin<2
    SVinview_least=10;
end

if size(sats,2)>1
    temp=find(sum(~isnan(sats)')>=SVinview_least);      % interested in the period with 10+ svs
else
    temp=find(~isnan(sats)');               % interested in the period of having data
end




tempt1=1;
tempt2=1;
t2=0;
t1=0;
for k=2:length(temp)
    if temp(k)==temp(k-1)+1
        tempt2=k;
    else
        if (tempt2-tempt1)>(t2-t1)
                t2=tempt2;
                t1=tempt1;
        end
            
        tempt1=k;
        tempt2=k;
        
    end
end

if (tempt2-tempt1)>(t2-t1)
        t2=tempt2;
        t1=tempt1;
end

t2=temp(t2);
t1=temp(t1);
