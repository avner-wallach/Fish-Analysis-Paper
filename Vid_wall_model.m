% make video showing tail effects on response to brass pole for example
% LFP, type 1 and type 2 cells
[ret,sname]=system('hostname');
if(strfind(sname,'Mac'))
    rootdir='/Volumes/aw3057/';
    sep='/';
else
    rootdir='Z:\';
    sep='\';
end
COL1=brewermap(10,'Set1');
COL=[226 139 138;210 149 0;74 133 34]/255;
nline=1;
wline=2;
bwidth=4;
msize=.5;
fontsize=7;
samplerate=3e4;
lp=100*2/samplerate;
[NN, Wn] = buttord( lp, lp .* [0.75], 3, 20);
[Bhp,Ahp] = butter(NN,Wn,'high');
%% load data

if(~exist('eod1'))
    load([rootdir,'analysis_data',sep,'20190131',sep,'output.mat']);
    eod1=eod;
    file1=file;
    ops1=ops;
end

load('Z:\mormyrid_data\fish_images\fish_silhouette.mat');
str_image=IMGS;
setenv('PXLSIZE','0.423'); %pxl to mm conversion
pxl2mm=0.423;
setenv('FRAMESCALE',num2str(ops.framescale));
setenv('SESSDATE','20190131');
setenv('DATAPATH','Z:\mormyrid_data');

lfpnum=4;
fdate='20190131';
% fdate='20210109';
dname=['Z:\mormyrid_data\',fdate,'\'];
fname=8; %file to take
% ss=29000;
ss=300;
framenum=300;
frameind=ss+[0:framenum];
% frameind=ss+[1:framenum];
if(~exist('fish_frames'))
    vidname=[dname,'video_',num2str(fname)];
    trackfile=[dname,'video_',num2str(fname),'_tracking'];
    datafile=[dname,'data_',num2str(fname)];    
    vid = VideoReader([vidname,'.avi']);
    load(datafile);
    [num,txt,raw] = xlsread([trackfile,'.csv']);
    F=figure;
    set(F,'Color',[0 0 0],'Position',[1 41 1920 964]);
    setenv('DATAPATH',ops.datapath);
    setenv('SESSDATE',fdate);
    setenv('FRAMESCALE',num2str(3));
    [fish_frames,ax_grid_f] = get_fish_image(fname,frameind,'headfix','off','track','off');
    [model_frames,ax_grid_m]=visualize_fish_modeling(vid,txt,num,ss,framenum,data.FILE.circle);

    X=data.FRAME.posture(frameind,1);
    Y=data.FRAME.posture(frameind,2);
    Az=data.FRAME.posture(frameind,3);
end
%% get wall position
rx=data.FILE.circle(3)/2; ry=data.FILE.circle(4)/2;    %ellipse radii
xc=data.FILE.circle(1) + rx; yc=data.FILE.circle(2) + ry; %center point
phi=atan2((Y-yc),(X-xc));    %azimuth in tank
R1=hypot((Y-xc),(X-yc)); %distance of fish from center
R2=(rx*ry)./sqrt((ry*cos(phi)).^2 + (rx*sin(phi)).^2); %distance of nearest point from cetner
%get wall allocentric position
xw=xc+R2.*cos(phi);  yw=yc+R2.*sin(phi);
%% get tail position
tail=nanmean(data.FRAME.posture(frameind,7:8),2);
%% get traces
cvar=13;    %number of channel to plot
ft=data.FRAME.t(frameind);
et=data.EOD.t;
clear frind Tr idx;
iii=find(et>=ft(1) & et<=ft(end));
e=1; frind=[]; idx=[];
while(e<=numel(iii))
    [m,i]=min(abs(ft-et(iii(e)))); %find frame closest to eod
    frind=[frind i];
    Tr(e,:)=data.EOD.traces(iii(e),1:180,cvar);
    e=e+1;
end
Trx=1+[1:size(Tr,2)]/30;
Tr=filtfilt(Bhp,Ahp,Tr')'/1e3;
%% generate fish images
K=51;
r = @(x) sqrt(x(:,1).^2 + x(:,2).^2);
w = @(x) atan2(x(:,2), x(:,1));
I0=imresize(str_image.cdata,20/17);
a=linspace(-1.35,1.35,K);
figure;
for i=1:numel(a)
    q = @(x) w(x) - a(i)*(x(:,2)>.0) .* (r(x).^2);
    f = @(x) [(r(x)) .* cos(q(x)), (r(x)) .* sin(q(x))];
    g = @(x, unused) f(x);

    tform2 = maketform('custom', 2, 2, [], g, []);
    I{i} = imtransform(I0, tform2, 'UData', [-1 1], 'VData', [-1 1]+.1, ...
        'XData', [-1 1], 'YData', [-1 1],'FillValues',255);
    imshow(I{i});
    drawnow;
end
%tail image for fig1
figure;
ind=[1 13 26 39 51]; 
q=linspace(1,0,numel(ind)+1);
for i=1:numel(ind)
    II=max(I{ind(i)},255*q(i+1));
    ai=image(II);
    set(ai,'AlphaData',(I{ind(i)}(:,:,1)<200));
    hold on;
end
axis('image');
set(gca,'Color','none','XColor','none','YColor','none');
%% get model images
[map_frames,Z,angle_v,R,TH,N,P1,P2]=get_tail_wall;
model_axes=get(gcf,'Children');
axis(model_axes(1),'image');
colormap(model_axes(1),invert_map(flipud(brewermap(64,'RdBu'))));
axis(model_axes(2),'image');
model_axes(2).CLim=[-0.75 0.75];
set(gcf,'Color',[0 0 0],'Units','normalized','Position',[0 0 1 1]);
map_slope=getframe(model_axes(1));
map_offset=getframe(model_axes(2));
%%
Coff=[-1 1];
Cslope=[-.75 .75];
bgcol=[0 0 0];
upsamp=2;
minsamp=3;
fontsize=14;

x_nbins=23;
y_nbins=23;

p=linspace(0,1,K+1);
% p=p(1:end-1);
p0=.2;
% p1=linspace(0,1-p0,K);
% p2=linspace(p0,1,K);
% p3=[p1;p2];
t=eod1(2).data(:,6);
r0=range(t)*p0;
% p=linspace(min(eod1(2).data(:,6)),max(eod1(2).data(:,6)),K+1);
q1=linspace(min(t),max(t)-r0,K);
q2=linspace(min(t)+r0,max(t),K);
q3=[q1;q2];

% q=quantile(unit1_eod.data(:,6),p);
%% LFP
% q=reshape(quantile(eod1(2).data(:,6),p3(:)),2,[]);
figure;
set(gcf,'Units','normalized','Position',[0 0 1 1],'Color',[0 0 0]);
for i=1:K
    plot_tuning_polar(eod1(2:4),file1(2),18,'t_col',6,'clim',[-2 2],'t_val',q3(:,i)',...
        'image',I{i},'fontsize',fontsize,'mfunc','mean','ops',ops1,'bgcol',bgcol,...
        'upsamp',.25,'r_nbins',25,'th_nbins',25,'ind_lim',[0e5 inf],'minsamp',15);
    colorbar('off');
%     set(gca,'Alim',[0 5]);
    frame_lfp(i)=getframe(gca);
end

plot_tuning_polar(eod1(2:4),file1(2),18,'mfunc','offset','t_col',6,'image',I0,...
        'clim',[-2 2],'bgcol',bgcol,...
        'upsamp',.5,'r_nbins',35,'th_nbins',55,'ind_lim',[0e5 inf],'minsamp',15);    
colorbar('off');
frame_lfp_offset=getframe(gca);

[Mlfp,Nlfp,hh]=plot_tuning_polar(eod1(2:4),file1(2),18,'mfunc','slope','t_col',6,'image',I0,...
        'clim',Cslope,'bgcol',bgcol,...
        'upsamp',.5,'r_nbins',25,'th_nbins',25,'ind_lim',[0e5 inf],'minsamp',15);
colorbar('off');
frame_lfp_slope=getframe(gca);
Cmap=get(gca,'Colormap');

%% compose video
frate=K/2;
d_2s=linspace(0,1,2*frate);
d_1s=linspace(0,1,frate);
d_05s=linspace(0,1,frate/2);
vidpath='C:\Users\sawte\Dropbox\efish_paper\Rev1\Videos\';
vidout=VideoWriter([vidpath,'WallachSawtell2023_MovieS2'],'mpeg-4');
vidout.FrameRate=frate;
vidout.Quality=25;
open(vidout);
Fall=figure;
set(Fall,'units','pixels','Position',[100 100 492*2 276*2],'Color',bgcol)
%% title page
title_image=imread('C:\Users\sawte\Dropbox\efish_paper\Videos\wall_vid_title.png');
A=axes('Position',[0 0 1 .95],'Color',[0 0 0]);
imshow(title_image);
A2=axes('Position',[0 0 1 .95],'Color','none','XColor','none','YColor','none');
R=rectangle('Position',[0 0 1 1],'LineStyle','none','FaceColor',[0 0 0 1]);
for i=1:frate
    R.FaceColor=[0 0 0 1-d_1s(i)];
    writeVideo(vidout,getframe(Fall));    
end
for i=1:frate*2
    writeVideo(vidout,getframe(Fall));    
end
for i=1:frate
    R.FaceColor=[0 0 0 d_1s(i)];
    writeVideo(vidout,getframe(Fall));    
end
delete(A);
delete(R);
%%
[ha, pos] = tight_subplot(2,4, [.01 .01],[.05 0.1],[.05 0],[],[.3 .3 .1 .3]); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
set(ha(:),'Color',bgcol,'Xcolor','none','Ycolor','none');
L0=annotation('textbox','String','0.5X Speed','Color',[1 1 1],'FontSize',fontsize,'Position',[0.01 0.4 0.2 0.1],'LineStyle','none');
%% fish movie
x0=403; y0=395; R=400;
[X,Y]=meshgrid(1:size(fish_frames(1).cdata,2),1:size(fish_frames(1).cdata,1));
blnkind=find(((X-x0).^2 + (Y-y0).^2)>R^2);
axes(ha(1));
cla;
j=1;
Im=fish_frames(j).cdata;
for d=1:3
    CD=Im(:,:,d);
    CD(blnkind)=0;
    Im(:,:,d)=CD;    
end
image(ax_grid_f{1},ax_grid_f{2},(Im));
axis('image');
set(gca,'XColor','none','YColor','none');
hold on;
H=plot(xw(j),yw(j),'o');
set(H,'MarkerFaceColor',COL1(1,:),'MarkerEdgeColor',[1 1 1],'MarkerSize',5);
H=plot([xc xw(j)],[yc yw(j)],':');
set(H,'LineWidth',2,'Color',[COL1(1,:) .25]);
Ocol=[1 1 1;0 0 0;COL1(1,:);[1 1 1]]
O(1)=annotation('textbox','String','Empty Tank','FontSize',fontsize,'Color',Ocol(1,:),'LineStyle','none',...
    'FontWeight','bold','HorizontalAlignment','center','Position',[0.11 0.92 0.15 0.05]);
O(2)=annotation('textarrow','String',['Implanted',char(10),'Fish'],'FontSize',fontsize,'Color',Ocol(2,:),...
    'HeadWidth',5,'HeadLength',5,'FontWeight','bold','HorizontalAlignment','center','Position',[.18 .73 -0.025 0]);
O(3)=annotation('textarrow','String','Wall position\newlineclosest to fish','FontSize',fontsize,'Color',Ocol(3,:),...
    'HeadWidth',5,'HeadLength',5,'FontWeight','bold','HorizontalAlignment','center','Position',[0.125 0.6411 -0.04 0.1094]);
O(4)=ylabel('Experiment','FontSize',fontsize*1.2,'Color',[1 1 1],'Position',[-150 600 1],...
    'FontWeight','bold','HorizontalAlignment','center','Visible','on');

A3=axes;
set(A3,'Position',pos{1}.*[1 1 .2 .25]+[pos{1}(3)*0 pos{1}(4)*.75 0 0],'Color','none','XColor',[1 1 1],'YColor',[1 1 1],'XAxisLocation','origin','YAxisLocation','left');
set(A3,'Xlim',[min(Trx) max(Trx)],'Ylim',[-1.5 max(Tr(:))],'YTick',[ -1 -.5 0],'Box','off','TickLength',[.05 .05],'TickDir','both','FontSize',fontsize*.5);
title('LFP','FontSize',fontsize,'Color',COL(1,:));

A2=axes('Position',[0 0 1 .95],'Color','none','XColor','none','YColor','none');
R=rectangle('Position',[0 0 1 1],'LineStyle','none','FaceColor',[0 0 0 1]);

hold on;
for i=1:frate    
    R.FaceColor=[0 0 0 1-d_1s(i)];
    for o=1:3
        O(o).Color=lighter(Ocol(o,:),d_1s(i)-1);
    end    
    writeVideo(vidout,getframe(Fall));    
end
for i=1:3*frate
    writeVideo(vidout,getframe(Fall));    
end
delete(O(1:3)); delete A2; delete R;
for j=1:numel(fish_frames)
    axes(ha(1));
    cla;
    Im=fish_frames(j).cdata;
    for d=1:3
        CD=Im(:,:,d);
        CD(blnkind)=0;
        Im(:,:,d)=CD;    
    end
    image(ax_grid_f{1},ax_grid_f{2},(Im));
    axis('image');
    set(gca,'XColor','none','YColor','none');
    hold on;
    H=plot(xw(j),yw(j),'o');
    set(H,'MarkerFaceColor',COL1(1,:),'MarkerEdgeColor',[1 1 1],'MarkerSize',5);
    ylabel('Experiment','FontSize',fontsize*1.2,'Color',[1 1 1],'Position',[-150 600 1],...
        'FontWeight','bold','HorizontalAlignment','center','Visible','on');
    %LFP traces
    axes(A3);    
    if(ismember(j,frind))
        cla(A3);
        trind=find(frind==j)+[-3:0];
        trind(trind<1)=[];
        H=plot(Trx,Tr(trind,:));
        set(H(1:end-1),'LineWidth',1,'Color',0.25*[1 1 1]);
        set(H(end),'LineWidth',2,'Color',COL(1,:));
        hold on;
        set(A3,'Color','none','XColor',[1 1 1],'YColor',[1 1 1],'XAxisLocation','origin','YAxisLocation','left');    
        set(A3,'Xlim',[min(Trx) max(Trx)],'Ylim',[-1.5 max(Tr(:))],'YTick',[ -1 -.5 0],'Box','off','TickLength',[.05 .05],'TickDir','both','FontSize',fontsize/2);
        title('LFP','FontSize',fontsize,'Color',COL(1,:));        
    end
    writeVideo(vidout,getframe(Fall));    
end    

%% LFP movie
A2=axes('Position',pos{2}.*[1 1 1 1.2],'Color','none','XColor','none','YColor','none');
R=rectangle('Position',[0 0 1 1],'LineStyle','none','FaceColor',[0 0 0 1]);
k=1;
for j=1:2
    for i=1:(2*K-1)
        axes(ha(2));
        imshow(frame_lfp(K-abs(i-K)).cdata);
        T=title({'Mean amplitude','(wall coordinates)'},'FontSize',fontsize,'FontName','Arial','FontWeight','bold','Color',[1 1 1]);        
        if(k<=numel(d_05s))
            axes(A2);
            R.FaceColor=[0 0 0 1-d_05s(k)];
            k=k+1;
        end        
        writeVideo(vidout,getframe(Fall));
    end
end

axes(ha(3));
ha(3).clo;
% imshow(frame_lfp_offset.cdata);
m=Mlfp(:);
n=Nlfp(:);
n(n<50)=nan;
n=n/max(n)*100+1;
S=scatter(1+randn(size(m))*.1,m,n,m,'filled','MarkerEdgeColor',0.5*[1 1 1],'MarkerFaceAlpha',1,'LineWidth',.1);
set(ha(3),'Xlim',[0 2],'Ylim',Cslope,'YTick',[-.5 0 .5 ],'Clim',Cslope,'XAxisLocation','origin','XTick',[]);
set(ha(3),'Colormap',Cmap,'Color','none','XColor',[1 1 1],'YColor',[1 1 1]);
T=title({'     Regression\newlineslope distribution'},...
    'FontSize',fontsize,'FontName','Arial','FontWeight','bold','Color',[1 1 1],...
    'HorizontalAlignment','center');    

axes(ha(4));
imshow(frame_lfp_slope.cdata);
T=title({'Tail reafference','map'},'FontSize',...
    14,'FontName','Arial','FontWeight','bold','Color',[1 1 1],...
    'HorizontalAlignment','center');        
A2=axes('Position',pos{3}.*[.945 1 8 1.3],'Color','none','XColor','none','YColor','none');
R=rectangle('Position',[0 0 1 1],'LineStyle','none','FaceColor',[0 0 0 1]);
k=1;

for i=1:2*frate
    if(k<=numel(d_05s))
        axes(A2);
        R.FaceColor=[0 0 0 1-d_05s(k)];
        k=k+1;
    end            
    writeVideo(vidout,getframe(Fall));
end 
%% model movie
axes(ha(5));
cla;
image(ax_grid_m{1},ax_grid_m{2},model_frames(1).cdata);
axis('image');
set(gca,'XColor','none','YColor','none');
O(1)=ylabel('Electrostatic Model','FontSize',fontsize*1.2,'Color',[1 1 1],...
    'FontWeight','bold','HorizontalAlignment','center','Visible','on','Position',[-150 600 1]);
A2=axes('Position',pos{5},'Color','none','XColor','none','YColor','none');
R=rectangle('Position',[0 0 1 1],'LineStyle','none','FaceColor',[0 0 0 1]);

for j=1:numel(model_frames)
    axes(ha(5));
    cla;
    image(ax_grid_m{1},ax_grid_m{2},model_frames(j).cdata);
    axis('image');
    set(gca,'XColor','none','YColor','none');
    ylabel('Electrostatic Model','FontSize',fontsize*1.2,'Color',[1 1 1],...
        'FontWeight','bold','HorizontalAlignment','center','Visible','on','Position',[-150 600 1]);
%     H=plot(xw(j),yw(j),'o');
%     set(H,'MarkerFaceColor',COL1(1,:),'MarkerEdgeColor',[1 1 1],'MarkerSize',5);
    writeVideo(vidout,getframe(Fall));    
end    
%% model results
A2=axes('Position',pos{6}.*[.945 1 1 1],'Color','none','XColor','none','YColor','none');
R=rectangle('Position',[0 0 1 1],'LineStyle','none','FaceColor',[0 0 0 1]);
k=1;
for j=1:2
    for i=1:(2*K-1)
        axes(ha(6));
        imshow(map_frames(K-abs(i-K)).cdata);
        if(k<=numel(d_05s))
            axes(A2);
            R.FaceColor=[0 0 0 1-d_05s(k)];
            k=k+1;
        end                    
        writeVideo(vidout,getframe(Fall));
    end
end


axes(ha(7));
ha(7).clo;
% imshow(map_offset.cdata);
m=P1(:);
n=N(:);
n=n/max(n)*25+1;
S=scatter(1+randn(size(m))*.1,m,n,m,'filled','MarkerEdgeColor',0.5*[1 1 1],'MarkerFaceAlpha',1,'LineWidth',.1);
set(ha(7),'Xlim',[0 2],'Ylim',Cslope,'YTick',[-.5 0 .5 ],'Clim',Cslope,'XAxisLocation','origin','XTick',[]);
set(ha(7),'Colormap',Cmap,'Color','none','XColor',[1 1 1],'YColor',[1 1 1]);

axes(ha(8));
imshow(map_slope.cdata);
A2=axes('Position',pos{7}.*[.945 1 8 1],'Color','none','XColor','none','YColor','none');
R=rectangle('Position',[0 0 1 1],'LineStyle','none','FaceColor',[0 0 0 1]);
k=1;

for i=1:2*frate;
    if(k<=numel(d_05s))
        axes(A2);
        R.FaceColor=[0 0 0 1-d_05s(k)];
        k=k+1;
    end            
    writeVideo(vidout,getframe(Fall)); 
end
close(vidout);