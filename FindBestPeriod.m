function [t1 t2]=FindBestPeriod(sats)

temp=sum(~isnan(sats)');

tempt1=1;
tempt2=1;
t2=0;
t1=0;
for k=2:length(temp)
    if k==length(temp) || temp(k)==temp(k-1)
        tempt2=k;
    else
        if temp(k-1)~=0 && (tempt2-tempt1)>(t2-t1)
                t2=tempt2;
                t1=tempt1;
        end
            
        tempt1=k;
        tempt2=k;
        
    end
end