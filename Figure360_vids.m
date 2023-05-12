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
bgcol=[1 1 1];
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
    set(F,'Color',bgcol,'Position',[1 41 1920 964]);
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
% [map_frames,Z,angle_v,R,TH,N,P1,P2]=get_tail_wall;
% model_axes=get(gcf,'Children');
% axis(model_axes(1),'image');
% colormap(model_axes(1),invert_map(flipud(brewermap(64,'RdBu'))));
% axis(model_axes(2),'image');
% model_axes(2).CLim=[-0.75 0.75];
% set(gcf,'Color',[0 0 0],'Units','normalized','Position',[0 0 1 1]);
% map_slope=getframe(model_axes(1));
% map_offset=getframe(model_axes(2));
%%
Coff=[-1 1];
Cslope=[-.75 .75];
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

%% compose video
frate=K/2;
d_2s=linspace(0,1,2*frate);
d_1s=linspace(0,1,frate);
d_05s=linspace(0,1,frate/2);
vidpath='C:\Users\sawte\Dropbox\efish_paper\Rev2 Final\videos\';
vidout=VideoWriter([vidpath,'Figure360_vid1'],'mpeg-4');
vidout.FrameRate=frate;
vidout.Quality=25;
open(vidout);
Fall=figure;
set(Fall,'units','pixels','Position',[100 100 276*2 276*2],'Color',bgcol)
%% fish movie
x0=403; y0=395; R=400;
[X,Y]=meshgrid(1:size(fish_frames(1).cdata,2),1:size(fish_frames(1).cdata,1));
blnkind=find(((X-x0).^2 + (Y-y0).^2)>R^2);
ha=axes();
pos{1}=ha.Position;
cla;
j=1;
Im=fish_frames(j).cdata;
for d=1:3
    CD=Im(:,:,d);
    CD(blnkind)=bgcol(1)*255;
    Im(:,:,d)=CD;    
end
image(ax_grid_f{1},ax_grid_f{2},(Im));
axis('image');
set(gca,'XColor','none','YColor','none');
hold on;
H=plot(xw(j),yw(j),'o');
set(H,'MarkerFaceColor',COL1(1,:),'MarkerEdgeColor',1-bgcol,'MarkerSize',15);
H=plot([xc xw(j)],[yc yw(j)],':');
set(H,'LineWidth',5,'Color',[COL1(1,:) .5]);
A3=axes;
set(A3,'Position',pos{1}.*[1 1 .2 .25]+[pos{1}(3)*0 pos{1}(4)*.75 0 0],'Color','none','XColor',1-bgcol,'YColor',1-bgcol,'XAxisLocation','origin','YAxisLocation','left');
set(A3,'Xlim',[min(Trx) max(Trx)],'Ylim',[-1.5 max(Tr(:))],'YTick',[ -1 -.5 0],'Box','off','TickLength',[.05 .05],'TickDir','both','FontSize',fontsize*.5);
writeVideo(vidout,getframe(Fall));    

for j=1:numel(fish_frames)
    axes(ha(1));
    cla;
    Im=fish_frames(j).cdata;
    for d=1:3
        CD=Im(:,:,d);
        CD(blnkind)=bgcol(1)*255;
        Im(:,:,d)=CD;    
    end
    image(ax_grid_f{1},ax_grid_f{2},(Im));
    axis('image');
    set(gca,'XColor','none','YColor','none');
    hold on;
    H=plot(xw(j),yw(j),'o');
    set(H,'MarkerFaceColor',COL1(1,:),'MarkerEdgeColor',1-bgcol,'MarkerSize',15);
    axes(A3);    
    if(ismember(j,frind))
        cla(A3);
        trind=find(frind==j)+[-3:0];
        trind(trind<1)=[];
        H=plot(Trx,Tr(trind,:));
        set(H(1:end-1),'LineWidth',1,'Color',0.75*[1 1 1]);
        set(H(end),'LineWidth',2,'Color',COL(1,:));
        hold on;
        set(A3,'Color','none','XColor',1-bgcol,'YColor',1-bgcol,'XAxisLocation','origin','YAxisLocation','left');    
        set(A3,'Xlim',[min(Trx) max(Trx)],'Ylim',[-1.5 max(Tr(:))],'YTick',[ -1 -.5 0],'Box','off','TickLength',[.05 .05],'TickDir','both','FontSize',fontsize/2);
    end
    writeVideo(vidout,getframe(Fall));    
end    
close(vidout);
%% model movie
vidout=VideoWriter([vidpath,'Figure360_vid2'],'mpeg-4');
vidout.FrameRate=frate;
vidout.Quality=100;
open(vidout);
axes(A3);    
delete(A3);
for j=1:numel(model_frames)
    axes(ha(1));
    cla;
    image(ax_grid_m{1},ax_grid_m{2},model_frames(j).cdata);
    axis('image');
    set(gca,'XColor','none','YColor','none');
    writeVideo(vidout,getframe(Fall));    
end
close(vidout);
%% traces
load('C:\Users\sawte\Dropbox\efish_paper\matlab\fish_frames.mat','p');
figure;
clim=[-1 1];
[Mz,N0,s1]=plot_tuning_polar(eod1(2),file1(2),18,'t_col',6,'clim',[-1 1],...
    'image',str_image.cdata,'fontsize',fontsize,'mfunc','slope','ops',ops1,'bgcol',bgcol,...
    'upsamp',0.25,'r_nbins',10,'th_nbins',15,'poly',p); 
x=s1.a_out(1).Children(1).XData;
J=[1:7];
Jd=linspace(1,7,K);
for i=1:3        
    for j=1:7
        y{i}(j,:)=s1.a_out(i).Children(j+(i==3)*2).YData;
    end
    yd{i}=flipud(interp1(J,y{i},Jd));
end

set(Fall,'Position',[100 100 286 200]);
ha(1).clo;
[hlfp,plfp]=tight_subplot(1,2,[.025 .05],[0 0],[0 0],[],[.5 .5],ha(1));
JJ=[1:K (K-1):-1:1];
for i=1:3
    vidout=VideoWriter([vidpath,'Figure360_vid4_',num2str(i)],'mpeg-4');
    vidout.FrameRate=frate;
    vidout.Quality=100;
    open(vidout);    
    for j=JJ
        axes(hlfp(1));    
        hlfp(1).clo;
        imshow(I{j});
        
        axes(hlfp(2));    
        hlfp(2).clo;
        H=plot(x,yd{i}(j,:));
        set(H,'Color',1-bgcol,'LineWidth',3);
        set(hlfp(2),'XLim',[0.0333 5],'YLim',[-2 1],'XColor','none','YColor','none','Color',bgcol);
        writeVideo(vidout,getframe(Fall));            
    end
    close(vidout);
end
%% generate field
plot_map('tank_radius',1e3,'wall_dist',5e2,'tail_p',0.6,'tail_angle',pi/7,'r_max',25,'bgcol',[1 1 1],'mpos','',...
    'plot_potential',0,'plot_field',1,'grid_M',100,'object_x',[10 3;-7 -15;6 -5],'object_R',[4 3 1],'object_c',[-1 -1 1])

plot_map('tank_radius',1e5,'wall_dist',5e4,'tail_p',0.6,'tail_angle',0,'r_max',35,'bgcol',[1 1 1],'mpos','','circular',0,...
    'plot_potential',0,'plot_field',1,'grid_M',50,'object_x',[10 5;-7 -30;8 -10],'object_R',[5 5 .75],'object_c',[-1 -1 1])