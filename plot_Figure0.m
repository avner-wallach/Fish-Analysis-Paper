% make panels for new Figure 1: Methodology, Afferent input in freely moving fish
%% methodology panels' data
adir='20190624'; seg=6; fnum=9; ss=13760; spkgroup=1; chid=1; rastind=[1:2];
% script for visualizing freely moving recordings
apath='Z:\analysis_data';
if(~exist('eod2'))
    load([apath,filesep,adir,filesep,'output'],'ops','eod','file');
    eod2=eod;
    ops2=ops;
    file2=file;
    load([apath,filesep,adir,filesep,adir,'_archieve'],'frame');
end
path(path,'Z:\GitHub\Fish-Model');
path(path,'Z:\GitHub\Fish-Pipeline')
if(~exist('eod1'))
    load('Z:\analysis_data\20190131\output.mat')
    eod1=eod;
    ops1=ops;
    file1=file;    
end
load('Z:\mormyrid_data\fish_images\fish_silhouette.mat');
str_image=IMGS;
setenv('FRAMESCALE',num2str(ops2.framescale));
setenv('SESSDATE','20190131');
setenv('DATAPATH','Z:\mormyrid_data');

setenv('PXLSIZE',num2str(ops2.pxlsize)); %pxl to mm conversion
pxl2mm=ops2.pxlsize;
fdate=num2str(ops2.seg(seg).dates(1));
dname=[ops2.datapath,'\',fdate,'\'];
fname=ops2.seg(seg).files(fnum); %file to take
cvar=ops2.seg(seg).LFPgroups{ops2.seg(seg).Spkgroups(spkgroup)};
datafile=[dname,'data_',num2str(fname)];
framenum=800; %number of frames
frameind=ss+[0:framenum];

%% get video and data files
if(~exist('data'))
    load(datafile);        %for traces
end

if(~exist('amp'))
    chns=find(ops2.outchans);
    ch=chns(cvar);
    sesspath=[ops2.datapath,'\',fdate,'\'];
    [t,amp]=read_bonsai_binary([sesspath,'amp_',num2str(fname)],ops2.samplerate,ops2.chan_num,ops2.blocksize,ch,'adc');
    t=t-data.FILE.offset;
    a=amp(:,1);
end
%% start figure
% bgc=[0 0 0];
bgc=[1 1 1];
fontsize=7;
nline=0.25;
wline=.5;
msize=3;

Fall=figure;
set(Fall,'Units','centimeters');
set(Fall,'OuterPosition', [2  2  11.4  10],'Color',bgc);

[ha, pos] = tight_subplot(2, 1, [.1 .0],[.075 .0],[.05 .05],[]); %Nh, Nw, [gap_h gap_w], [lower upper], [left right] 
[h1,p1]=tight_subplot(1, 2, [.0 .0],[.0 .00],[.00 .0],[],[.7 .3],ha(1)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
[h2,p2]=tight_subplot(1, 2, [.0 .05],[.0 .0],[.0 .0],[],[.5 .5],ha(2)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
set(h1([1]),'Visible','off'); %schemes
%%
COL=colormap(ha(1),'lines');
% COL(4,:)=[236 157 118]/255; %afferent color
COL(4,:)=[226 139 138]/255; %afferent color
COL(5,:)=[64 128 0]/255;%[57 181 74]/255; output color
COL(6,:)=[210 148 37]/255;%COL(7,:);%interneuron
COL(7,:)=[86 124 141]/255;%output color
CHSL=rgb2hsl(COL);

ft=data.FRAME.t(frameind);
acc(:,1)=data.FRAME.accx;
acc(:,2)=data.FRAME.accy;
acc(:,3)=data.FRAME.accz;
[coef,score,latent,tqs,explained]=pca(acc);
posture=data.FRAME.posture;
rfin=posture(:,9)/2; rfin=rfin-nanmedian(rfin);
lfin=posture(:,10)/2; lfin=lfin-nanmedian(lfin);
tpos=mean(posture(:,7:8),2);
x1=posture(:,1); x1=(x1-nanmin(x1))/(nanmax(x1)-nanmin(x1));
y1=posture(:,2); y1=(y1-nanmin(y1))/(nanmax(y1)-nanmin(y1));
az=posture(:,3)/pi;
eind=[eod2(seg).ind0(fnum):eod2(seg).ind0(fnum+1)-1];   %eod indices for file
edata=eod2(seg).data(eind,:);   %eod data for file
et=eod2(seg).t(eind,:);    %eod2 time

%% plot ephys data
xlim_=ft(1)+[13.71 13.86];
ind=find(inrange(t,xlim_));
a1=a(ind);
axes(h2(1));
h2(1).clo;
H=plot(t(ind)-xlim_(1),a1,'Color',0.5*[1 1 1],'LineWidth',nline/2);
hold on;
%artifact
ett=et(inrange(et,xlim_));
ind=find(ismember(t,ett))'; %eod indices in t
Ix=[-15:1:30];
Idx=ind*ones(size(Ix))+ones(size(ind))*Ix;    
H=plot(t(Idx')-xlim_(1),a(Idx'),'Color',1-bgc,'LineWidth',wline);
%LFP
Ix=[30:1:6*30];
Idx=ind*ones(size(Ix))+ones(size(ind))*Ix;    
H=plot(t(Idx')-xlim_(1),a(Idx'),'Color',COL(4,:),'LineWidth',wline);
%spikes
Ix=[-15:1:45];
Rt=double(data.RAST(:,1))/ops2.samplerate-data.FILE.offset;
ind73=find(inrange(Rt,ft([1 end])) & data.RAST(:,2)==73);
ind=find(ismember(t,Rt(ind73)))'; %eod indices in t
Idx=ind*ones(size(Ix))+ones(size(ind))*Ix;    
Tr73=a(Idx);
H=plot(t(Idx')-xlim_(1),Tr73','Color',COL(6,:),'LineWidth',wline);
ind72=find(inrange(Rt,ft([1 end])) & data.RAST(:,2)==72);
ind=find(ismember(t,Rt(ind72)))'; %eod indices in t
Idx=ind*ones(size(Ix))+ones(size(ind))*Ix; 
Tr72=a(Idx);
H=plot(t(Idx')-xlim_(1),Tr72','Color',COL(5,:),'LineWidth',wline);
set(gca,'Xlim',[0 0.15],'Ylim',[-500 1000],'FontSize',fontsize);
set(gca,'XTick',[0 0.05 0.1 0.15],'XTickLabel',[0 50 100 150],...
    'YTick',[-500 0 500 1e3],'YTickLabel',{'-0.5','0','0.5','1'},'Box','off')
set(h2(1),'TickLength',[0.02 0.02],'TickDir','both');

% insets
A6=axes();
A6.clo;
dx=1;
set(A6,'Position',h2(1).Position.*[1 1 .55 .5] + [h2(1).Position(3)*.25 h2(1).Position(4)*.5 0 0]);
x=[1:size(eod2(seg).avtraces,2)]/ops2.samplerate*1e3+ops2.traceblnk;
ind=find(inrange(x,[1 6]));
Y=eod2(seg).avtraces(:,:,ch(1),1);
plot_graded(x(ind),Y(2:end,ind),flipud(colormap_graded(COL(4,:),0.375,1,size(Y,1)-1)));
H1=get(A6,'Children');
set(H1,'LineWidth',wline);
hold on;
H=plot([4 5],[200 200],'k',[4 4],[200 300],'k');
set(H,'LineWidth',wline);
x0=x(ind(end))+dx;
x=Ix/30+x0;
Tr73=Tr73-median(Tr73(:,1:10),2)*ones(1,size(Tr73,2));
H=plot(x,Tr73','Color',lighter(COL(6,:),0.85),'LineWidth',nline);
hold on;
plot(x,median(Tr73,1),'Color',COL(6,:),'LineWidth',wline);
Tr72=Tr72-median(Tr72(:,1:10),2)*ones(1,size(Tr72,2));
x0=x(end)+dx;
x=Ix/30+x0;
H=plot(x,Tr72','Color',lighter(COL(5,:),0.85),'LineWidth',nline);
hold on;
plot(x,median(Tr72,1),'Color',COL(5,:),'LineWidth',wline);
set(A6,'XTick',[],'YTick',[],'Xlim',[1-dx/2 x(end)+dx/2],'Ylim',[-550 750],'Box','on');

%% plot behavioral variables
axes(h1(2));
h1(2).clo;
[hb,pb]=tight_subplot(3, 1, [.0 .0],[.0 .00],[.00 .0],[.2 .5 .3],[],h1(2)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
axes(hb(3));
hb(3).clo;
k=0;
H=plot(ft-ft(1),x1(frameind)+k);
hold on;
set(H,'LineWidth',wline,'Color',COL(7,:));
k=k+1;
H=plot(ft-ft(1),y1(frameind)+k);
set(H,'LineWidth',wline,'Color',COL(7,:));
k=k+1;
H=plot(ft-ft(1),az(frameind)+k);
set(H,'LineWidth',wline,'Color',COL(7,:));
xlim=[0 ft(end)-ft(1)];
set(hb(3),'FontSize',fontsize,'YTick',[],'Xlim',xlim,'XTick',[0 5 10 15],'Ylim',[0 3.25],'Box','off');
set(hb(3),'TickLength',[0.025 0],'TickDir','both');

axes(hb(2));
hb(2).clo;
k=0;
H=plot(ft-ft(1),lfin(frameind)+k);
hold on;
set(H,'LineWidth',wline,'Color',COL(3,:));
k=k+1;
H=plot(ft-ft(1),rfin(frameind)+k);
set(H,'LineWidth',wline,'Color',COL(3,:));
k=k+1;
H=plot(ft-ft(1),score(frameind,2)+k);
set(H,'LineWidth',wline,'Color',COL(11,:));
k=k+1;
H=plot(ft-ft(1),score(frameind,1)+k);
set(H,'LineWidth',wline,'Color',COL(11,:));
k=k+1;
H=plot(ft-ft(1),tpos(frameind)+k);
set(H,'LineWidth',wline,'Color',COL(1,:));
xlim=[0 ft(end)-ft(1)];
set(hb(2),'FontSize',fontsize,'YTick',[],'Xlim',xlim,'Ylim',[-.5 4.5],'XTick',[],'Box','off');

% [hc,pc]=split_axes(hb(1),4,1,0.0,0.005);
isi=diff(data.EOD.t);
et=data.EOD.t(inrange(data.EOD.t,ft([1 end])));
ifr=1./isi(inrange(data.EOD.t,ft([1 end])));
axes(hb(1));
H=plot(et-ft(1),ifr,'.');
set(H,'MarkerSize',msize,'Color',COL(2,:),'MarkerFaceColor',COL(2,:));
set(hb(1),'FontSize',fontsize,'Xlim',xlim,'Box','off','XTick',[],'Ylim',[4 12],'YTick',[5 10],'TickLength',[0.02 0.02],'TickDir','both');
%% rasters
% lfp_eod=eod1(4); lfp_file=file1(4); lfp_col=17; %lfp_zscale=[430 40];
lfp_eod=eod2(6); lfp_file=file2(6); lfp_col=21;
unit1_eod=eod2(6); unit1_file=file2(6); unit1_col=24; unit1_ind=ops2.seg(6).ind_lim(2,:); unit1_clim=[0 2]; unit1_lfpcol=21; unit1_rast=unit1_eod.raster{2};lfp_chan=13;
% unit2_eod=eod1(4); unit2_file=file1(4); unit2_col=20; unit2_ind=ops1.seg(4).ind_lim(2,:); unit2_clim=[0 2]; unit2_lfpcol=17; unit2_rast=unit2_eod.raster{2};
unit2_eod=eod2(6); unit2_file=file2(6); unit2_col=23; unit2_ind=ops2.seg(6).ind_lim(1,:); unit2_clim=[0 1.5]; unit2_lfpcol=21; unit2_rast=unit2_eod.raster{1};

[h21,pos21] = tight_subplot(1,2,[0 0],[0 0],[0 .0],[],[.25 .75],h2(2));
% h1(1).Visible='off';
[hL, posL] = tight_subplot(2, 1, [.0 .0],[.0 0.0],[.0 0.0],[.5 .5],[],h21(1)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
[hp, posp] = tight_subplot(2, 2, [.0 .03],[.0 0.0],[.0 0],[.5 .5],[],h21(2)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]

%LFP traces
avtraces=unit1_eod.avtraces(2:2:end,:,lfp_chan,1)/500;

cmp=flipud(colormap_graded(COL(4,:),.375,.85,size(avtraces,1)));
if(bgc==[0 0 0])
    cmp=flipud(cmp);
end
x=ops1.traceblnk+[1:size(unit1_eod.avtraces,2)]/30;
axes(hL(1));
hL(1).clo;
for i=1:size(avtraces,1)
    H=plot(x,avtraces(i,:)+i);
    set(H,'Color',cmp(i,:),'LineWidth',wline);
    hold on;        
end
x0=3.5; y0=i+.25;
H=plot(x0+[0 1],y0*[1 1],'k',x0*[1 1],y0+[0 1],'k');
set(H,'LineWidth',nline);
set(hL(1),'Xlim',[1 6],'Ylim',[0.5 i+1.5],'YColor','none','XColor','none','Color','none');
hL(2).Visible='off';

cmap{1}=hsl2rgb([ones(64,1)*CHSL(6,1:2) linspace(.95,0.25,64)']);
if(bgc==[0 0 0])
    cmap{1}=flipud(cmap{1});
end
[hr,hs]=plot_sorted_raster(unit1_rast,unit1_eod.data(:,unit1_lfpcol),ops2,unit1_ind,5e2,5,.5,cmap,bgc);
% close(gcf);
set(hs(2),'YLim',[0 200]);
set(hs(2),'Xtick',[-20 0 20 40],'YTick',[50 100 150],'TickLength',[0.02 .02],'TickDir','Both','FontSize',fontsize,'Color','none');
set(hs(1),'XColor','none','Color','none');
for i=1:numel(hs)
    cb=get(hs(i),'Children');
    for j=1:numel(cb)
        if(strcmp(cb(j).Type,'line'))
            cb(j).LineWidth=wline;
        elseif(strcmp(cb(j).Type,'Scatter'))
            cb(j).SizeData=msize/2;
        end
    end
    set(hs(i),'XLim',[-10 40]);    
    copyaxes(hs(i),hp(2*i-1),1);
    hp(2*i-1).Colormap=hs(i).Colormap;
end
set(hp(3),'TickLength',[0.05 .02],'TickDir','Both','FontSize',fontsize,'Color','none');
cmap{1}=hsl2rgb([ones(64,1)*CHSL(5,1:2) linspace(.95,0.25,64)']);
if(bgc==[0 0 0])
    cmap{1}=flipud(cmap{1});
end
[hr,hs]=plot_sorted_raster(unit2_rast,unit2_eod.data(:,unit2_lfpcol),ops2,unit2_ind,5e2,5,.5,cmap,bgc);
% close(gcf);
set(hs(2),'YLim',[0 80]);
set(hs(2),'Xtick',[-20 0 20 40],'YTick',[20 40 60],'TickLength',[0.02 .02],'TickDir','Both','FontSize',fontsize,'Color','none');
set(hs(1),'XColor','none','Color','none');

for i=1:numel(hs)
    cb=get(hs(i),'Children');
    for j=1:numel(cb)
        if(strcmp(cb(j).Type,'line'))
            cb(j).LineWidth=wline;
        elseif(strcmp(cb(j).Type,'scatter'))
            cb(j).SizeData=msize/2;
        end
    end    
    set(hs(i),'XLim',[-10 40]);
    copyaxes(hs(i),hp(i*2),1);
    hp(2*i).Colormap=hs(i).Colormap;
end
set(hp(3),'TickLength',[0.05 .02],'TickDir','Both','FontSize',fontsize,'Color','none');
aaa=findobj('Type','Area');
set(aaa(:),'ShowBaseline','off');

% axes(hp(4));
% x0=25; y0=35;
% H=plot(x0+[0 10],y0*[1 1],'k',x0*[1 1],y0+[0 40],'k');
% set(H,'LineWidth',nline);
% N=2e4;
% Ks=50;
% axes(h2(3));
% h2(3).clo;
% % h2(3).Position=pos{9}+[0.035 0.04 -0.025 -0.045];
% ind=randperm(size(unit1_eod.data,1),N);
% x=unit1_eod.data(ind,unit1_lfpcol);
% y=(unit1_eod.data(ind,unit1_col));
% y=nanzscore(y);
% y1=smooth(x,y,Ks,'moving')+.1*randn(size(y));
% ind=find(~isnan(x+y));
% f1=fit(x(ind),y(ind),'poly1');
% 
% X=sort(x(~isnan(x)));
% p = predint(f1,X,0.95,'functional','on');
% % A1=area(x,p(:,2),-1,'LineStyle','none','FaceColor',COL(2,:),'FaceAlpha',0.25);
% A1=patch([X;flipud(X)],[p(:,2);flipud(p(:,1))],COL(6,:),'LineStyle','none','FaceColor',COL(5,:),'FaceAlpha',0.25);
% hold on;
% H=plot(f1); set(H,'Color',COL(6,:),'LineWidth',wline); legend('off');
% % A2=area(x,p(:,1),-1,'LineStyle','none','FaceColor',[1 1 1],'FaceAlpha',1);
% scatter(x,y1,msize,'filled','MarkerFaceAlpha',0.05,'MarkerFaceColor',COL(6,:));
% 
% ind=randperm(size(unit2_eod.data,1),N);
% x=unit2_eod.data(ind,unit2_lfpcol);
% y=(unit2_eod.data(ind,unit2_col));
% y=nanzscore(y);
% y1=smooth(x,y,Ks,'moving')+.1*randn(size(y));
% ind=find(~isnan(x+y));
% f2=fit(x(ind),y(ind),'poly1');
% X=sort(x(~isnan(x)));
% p = predint(f2,X,0.95,'functional','on');
% A1=patch([X;flipud(X)],[p(:,2);flipud(p(:,1))],COL(5,:),'LineStyle','none','FaceColor',COL(5,:),'FaceAlpha',0.25);
% % A1=area(x,p(:,2),-1,'LineStyle','none','FaceColor',COL(3,:),'FaceAlpha',0.25);
% hold on;
% H=plot(f2); set(H,'Color',COL(5,:),'LineWidth',wline); legend('off');
% % A2=area(x,p(:,1),-1,'LineStyle','none','FaceColor',[1 1 1],'FaceAlpha',1);
% scatter(x,y1,msize,'filled','MarkerFaceAlpha',0.05,'MarkerFaceColor',COL(5,:),'MarkerEdgeAlpha',0);
% set(gca,'TickLength',[.025 .025],'TickDir','both');
% set(gca,'Ylim',[-.75 1],'Xlim',[-1 2.5],'FontSize',fontsize,'box','off','XLabel',[],'YLabel',[],...
%     'Xtick',[-1 0 1 2],'YTick',[-.5 0 .5 1],'XTickLabelMode','auto','YTickLabelMode','auto');

%% remove images before copy
% hs={h2(2),h31(1),h31(2),h33(1),h33(2)};
srf=findobj('Type','Surface');
ptch=findobj('Type','Patch');
img=findobj('Type','Image');
set(srf,'Visible','off');
set(img,'Visible','off');
for i=1:numel(ptch)
    if(ptch(i).Parent==h3(5))
       set(ptch(i),'Visible','off');
    end
end
print('-clipboard','-dmeta');
%% remove all but figures before copy to bmp
hall=findobj();
for i=1:numel(hall)
    if(~strcmp(hall(i).Type,'root') & ~strcmp(hall(i).Type,'figure'))
        set(hall(i),'Visible','off');
    end
end
srf=findobj('Type','Surface');
ptch=findobj('Type','Patch');
img=findobj('Type','Image');
set(srf,'Visible','on');
set(img,'Visible','on');
for i=1:numel(ptch)
    if(ptch(i).Parent==h3(5))
       set(ptch(i),'Visible','on');
    end
end
print('-clipboard','-dbitmap','-r720'); %copy maps

        
