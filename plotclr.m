function plotclr(x,y,v,prn,Plot_ValueRange,marker)
% clement.ogaja@gmail.com
    global Sat_Capacity     % 1-32 : GPS ; 33-64 : GLONASS ; 65-96: GALILEO
    global SatList

if nargin <6
    marker='.';
end

% map=colormap;
temp=hsv(Sat_Capacity*3);temp=temp(1:Sat_Capacity,:);
map=colormap(temp);   % create colormap as red-->yellow-->green
% map=colormap(hsv(Sat_Capacity*6));map=map(1:Sat_Capacity,:);    % create colormap as red-->yellow-->green

if Plot_ValueRange(1)>Plot_ValueRange(2)
    map=flipud(map);
end


[row col]=size(x);

% Plot the points
tempshift=0;
if Plot_ValueRange(1)<0
    v=v-Plot_ValueRange(1);
    
    tempshift=Plot_ValueRange(1);
    Plot_ValueRange=Plot_ValueRange-Plot_ValueRange(1);
end

    miv=min(Plot_ValueRange);
    mav=max(Plot_ValueRange);


hold on
for j=1:col,
    for i=1:row%length(x)
        in=round((v(i,j)-miv)*(length(map)-1)/(mav-miv));
        %--- Catch the out-of-range numbers
        if in<=0;in=1;end
        if in > length(map);in=length(map);end
        if isnan(in);in=1;end
%         [i j ]
%         v(i,j)
%         map(in,:)
        if v(i,j)>0
            plot3(x(i,j),y(i,j),v(i,j),marker,'color',map(in,:),'MarkerFaceColor',map(in,:))
        end
        if i==row,
            ht=text(x(i,j),y(i,j),SatList{prn(j)'});set(ht,'fontsize',8);            
        end
    end
end
hold off

% Re-format the colorbar
h=colorbar;
colormark_num=5;
set(h,'fontsize',8);
set(h,'ylim',[0 1])
yal=linspace(0,1,colormark_num);
set(h,'ytick',yal);
% Create the yticklabels
ytl=linspace(miv,mav,colormark_num)+tempshift;
s=char(colormark_num,4);
for i=1:colormark_num
    if min(abs(ytl)) >= 0.001
        B=sprintf('%4.2f',ytl(i));
    else
        B=sprintf('%4.2E',ytl(i));
    end
    s(i,1:length(B))=B;
end
set(h,'yticklabel',s);
% 
% set(h,'fontsize',8);
% % set(get(h,'ylabel'),'string','[m]');
% set(h,'ylim',[1 length(map)]);
% yal=linspace(1,length(map),6);
% set(h,'ytick',yal);
% % Create the yticklabels
% ytl=linspace(miv,mav,6)+tempshift;
% s=char(6,4);
% for i=1:6
%     if min(abs(ytl)) >= 0.001
%         B=sprintf('%4.2f',ytl(i));
%     else
%         B=sprintf('%4.2E',ytl(i));
%     end
%     s(i,1:length(B))=B;
% end
% set(h,'yticklabel',s);

grid on
view(2)

