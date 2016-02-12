function plotclr(x,y,v,prn,color_define,marker)
% clement.ogaja@gmail.com

if nargin <6
    marker='.';
end

map=colormap;


[row col]=size(x);

% Plot the points
tempshift=0;
if color_define(1)<0
    v=v-color_define(1);
    
    tempshift=color_define(1);
    color_define=color_define-color_define(1);
end

    miv=color_define(1);
    mav=color_define(2);


hold on
for j=1:col,
    for i=1:row%length(x)
        in=round((v(i,j)-miv)*(length(map)-1)/(mav-miv));
        %--- Catch the out-of-range numbers
        if in==0;in=1;end
        if in > length(map);in=length(map);end
        if isnan(in);in=1;end
%         [i j ]
%         v(i,j)
%         map(in,:)
        if v(i,j)>0
            plot3(x(i,j),y(i,j),v(i,j),marker,'color',map(in,:),'markerfacecolor',map(in,:))
        end
        if i==row,
            if prn(j) < 10,
                ht=text(x(i,j),y(i,j),[' S0',num2str(prn(j))]);set(ht,'fontsize',8);
            else
                ht=text(x(i,j),y(i,j),[' S',num2str(prn(j))]);set(ht,'fontsize',8);
            end
        end
    end
end
hold off

% Re-format the colorbar
h=colorbar;

set(h,'fontsize',8);
set(get(h,'ylabel'),'string','[m]');
set(h,'ylim',[1 length(map)]);
yal=linspace(1,length(map),6);
set(h,'ytick',yal);
% Create the yticklabels
ytl=linspace(miv,mav,6)+tempshift;
s=char(6,4);
for i=1:6
    if min(abs(ytl)) >= 0.001
        B=sprintf('%4.2f',ytl(i));
    else
        B=sprintf('%4.2E',ytl(i));
    end
    s(i,1:length(B))=B;
end
set(h,'yticklabel',s);
grid on
view(2)

