function bandview(theta,rho,z,prn,Plot_ValueRange)

%BANDVIEW  sky image plot in a bandview.
%   BANDVIEW(THETA, RHO, Z, PRN,Plot_ValueRange) makes a geographically oriented plot  
%   in a band view of the angle THETA, in radians, versus the radius  RHO.
%
%   Lei Yang

if isstr(theta) | isstr(rho)
    error('Input arguments must be numeric.');
end
if ~isequal(size(theta),size(rho))
    error('THETA and RHO must be the same size.');
end

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   get(cax, 'FontSize'), ...
    'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
    'DefaultTextUnits','data')

% only do grids if hold is off
if ~hold_state

% make a radial grid
    hold on;
    maxrho = max(abs(rho(:)));
    hhh=plot([-maxrho -maxrho maxrho maxrho],[-maxrho maxrho maxrho -maxrho]);
%     set(gca,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
    v = [get(cax,'xlim') get(cax,'ylim')];
%     v=[-1 1 -1 1];
    ticks = sum(get(cax,'ytick')>=0);
    delete(hhh);
% check radial limits and ticks
    rmin = 0; rmax = v(4); rticks = max(ticks-1,2);
    if rticks > 5   % see if we can reduce the number
        if rem(rticks,2) == 0
            rticks = rticks/2;
        elseif rem(rticks,3) == 0
            rticks = rticks/3;
        end
    end

% define a circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = 1:(length(th)-1)/4:length(th);
    xunit(inds(2:2:4)) = zeros(2,1);
    yunit(inds(1:2:5)) = zeros(3,1);
% plot background if necessary
    if ~isstr(get(cax,'color')),
       patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
             'edgecolor',tc,'facecolor',get(gca,'color'),...
             'handlevisibility','off');
    end

% draw radial circles
    c82 = cos(82*pi/180);
    s82 = sin(82*pi/180);
    rinc = (rmax-rmin)/rticks;
    for i=(rmin+rinc):rinc:rmax
        hhh = plot(xunit*i,yunit*i,ls,'color',tc,'linewidth',1,...
                   'handlevisibility','off');
        %text((i+rinc/20)*c82,(i+rinc/20)*s82, ...
        %    ['  ' num2str(i)],'verticalalignment','bottom',...
        %    'handlevisibility','off')
    end
    set(hhh,'linestyle','-') % Make outer circle solid

% plot spokes
    th = (1:6)*2*pi/12;
    cst = sin(th); snt = cos(th);
    cs = [-cst; cst];
    sn = [-snt; snt];
    plot(rmax*cs,rmax*sn,ls,'color',tc,'linewidth',1,...
         'handlevisibility','off')

% annotate spokes in degrees
    rt = 1.1*rmax;
    for i = 1:length(th)
        text(rt*cst(i),rt*snt(i),int2str(i*30),...
             'horizontalalignment','center',...
             'handlevisibility','off');
        if i == length(th)
            loc = int2str(0);
        else
            loc = int2str(180+i*30);
        end
        text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center',...
             'handlevisibility','off')
    end

% set view to 2-D
    view(2);
% set axis limits
    axis(rmax*[-1 1 -1.15 1.15]);
end

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits',fUnits );

% transform data to Cartesian coordinates.
xx = rho.*sin(theta);
yy = rho.*cos(theta);

% plot data on top of grid
if nargin <5
    Plot_ValueRange=[min(min(z)) max(max(z))];
end
plotclr(xx,yy,z,prn,Plot_ValueRange);
        
if ~hold_state
    set(gca,'dataaspectratio',[1 1 1]), axis off; set(cax,'NextPlot',next);
end
set(get(gca,'xlabel'),'visible','on')
set(get(gca,'ylabel'),'visible','on')
end

function plotsquare(x,y,v,prn,Plot_ValueRange,marker)

x=THETA/pi*180;
y=(1-RHO)*90;
v=sample_SatVal_Active;
 prn=SV_visiable;

    global Sat_Capacity     % 1-32 : GPS ; 33-64 : GLONASS ; 65-96: GALILEO
    global SatList

    if nargin <6,    marker='.';    end

    % create colormap as red-->yellow-->green
    temp=hsv(32*3);temp=temp(1:32,:);
    map=colormap(temp);   

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

            if v(i,j)>0
                plot(x(i,j),y(i,j),marker,'color',map(in,:),'markerfacecolor',map(in,:))
            end
            if i==row,
                ht=text(x(i,j),y(i,j),SatList{prn(j))'};set(ht,'fontsize',8);
            end
        end
    end
    hold off
    
    set(gca,'xlabel','Azimuth (degree)')
    set(gca,'ylabel','Elevation (degree)')

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

grid on
end
