% make panels for new Figure 1: Methodology, Afferent input in freely moving fish
[ret,sname]=system('hostname');
if(strfind(sname,'Mac'))
    rootdir='/Volumes/aw3057/';
    sep='/';
else
    rootdir='Z:\';
    sep='\';
end
%% methodology panels' data
adir='20190624'; seg=6; fnum=9; ss=13760; spkgroup=1; chid=1; rastind=[1:2];
% script for visualizing freely moving recordings
apath='Z:\analysis_data';
if(~exist('eod1'))
    load([apath,filesep,adir,filesep,'output'],'ops','eod','file');
    eod1=eod;
    ops1=ops;
    file2=file;
    load([apath,filesep,adir,filesep,adir,'_archieve'],'frame');
end
if(~exist('eod2'))
    load([rootdir,'analysis_data',sep,'20200113',sep,'output.mat']);
    eod2=eod;
    file2=file;
    ops2=ops;
end

path(path,'Z:\GitHub\Fish-Model');
path(path,'Z:\GitHub\Fish-Pipeline')
load('Z:\mormyrid_data\fish_images\fish_silhouette.mat');
str_image=IMGS;
setenv('FRAMESCALE',num2str(ops1.framescale));
setenv('SESSDATE','20190131');
setenv('DATAPATH','Z:\mormyrid_data');

setenv('PXLSIZE',num2str(ops1.pxlsize)); %pxl to mm conversion
pxl2mm=ops1.pxlsize;
fdate=num2str(ops1.seg(seg).dates(1));
dname=[ops1.datapath,'\',fdate,'\'];
fname=ops1.seg(seg).files(fnum); %file to take
cvar=ops1.seg(seg).LFPgroups{ops1.seg(seg).Spkgroups(spkgroup)};
datafile=[dname,'data_',num2str(fname)];
framenum=800; %number of frames
frameind=ss+[0:framenum];

%% get video and data files
if(~exist('data1'))
    load(datafile);        %for traces
    data1=data;
end

if(~exist('amp'))
    chns=find(ops1.outchans);
    ch=chns(cvar);
    sesspath=[ops1.datapath,'\',fdate,'\'];
    [t,amp]=read_bonsai_binary([sesspath,'amp_',num2str(fname)],ops1.samplerate,ops1.chan_num,ops1.blocksize,ch,'adc');
    t=t-data1.FILE.offset;
    a=amp(:,1);
end
%% rec dynamics 
Q=[.1 .25 .5 .75 .9];
K=100;
lfpind=[7,8];
lfpcol=18; scaleind=[4 4 5];
unit1_name=[32 33 33];
% unit2_name=[0 32 32];
unit2_name=[42 42 42];

unit1col=[19 21 21];
% unit2col=[0 20 20];
unit2col=[20 22 22];
% segs=[3:5];
segs=[3:4];
Tr=[];
U1=[];
U2=[];
dt=[];
dtswitch=[];
ISI=[];
for i=1:numel(segs)
    clear lfp_q unit1_q;
    for j=1:numel(eod2(segs(i)).ind0)-1
        ind=[eod2(segs(i)).ind0(j):eod2(segs(i)).ind0(j+1)-1];
        lfp=eod2(segs(i)).data(ind,lfpcol)*eod2(segs(i)).LFPs(scaleind(i))+eod2(segs(i)).LFPm(scaleind(i));
%         lfp_q(j,:)=nanmean(lfp)+nanstd(lfp)*[-1 0 1];
        lfp_q(j,:)=quantile(lfp,Q);
%         unit1_q(j,:)=nanmean(unit1)+nanstd(unit1)*[-1 0 1];;        
    end
    isi=eod2(segs(i)).data(:,1);
    isi(isnan(isi))=nanmean(isi);        
    if(unit1col(i))
        unit1=smooth(eod2(segs(i)).data(:,unit1col(i)),K);
    else
        unit1=zeros(size(isi));
    end
    if(unit2col(i))
        unit2=smooth(eod2(segs(i)).data(:,unit2col(i)),K);
    else
        unit2=zeros(size(isi));
    end
%     Tr=cat(3,Tr,permute(eod2(i).avtraces(:,:,lfpind,:),[1 2 4 3]));
    Tr=[Tr;lfp_q];
    U1=[U1;unit1];
    U2=[U2;unit2];
    ISI=[ISI;isi];
    dt=[dt;eod2(segs(i)).t(eod2(segs(i)).ind0(2:end)-1)/3600];
    dtswitch=[dtswitch;sum(eod2(segs(i)).t(eod2(segs(i)).ind0(2:end)-1)/3600)];
end
T=cumsum(dt);
Te=cumsum(ISI)/3600;
tswitch=[1;cumsum(dtswitch)];

%% start figure
% bgc=[0 0 0];
bgc=[1 1 1];
fontsize=7;
nline=0.25;
wline=.5;
msize=3;
tlength=0.0125;

Fall=figure;
set(Fall,'Units','centimeters');
set(Fall,'OuterPosition', [2  2  12  14],'Color','none');

[ha, pos] = tight_subplot(3, 1, [.075 .0],[.05 .0],[.05 .0],[.33 .34 .33]); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
% first row
[h1,p1]=tight_subplot(1, 2, [.0 .075],[.0 .00],[.00 .0],[],[.66 .34],ha(1)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
% second row
[h2,p2]=tight_subplot(1, 3, [.0 .075],[.0 .00],[.00 .0],[],[.2 .3 .5],ha(2)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]

set(h1([1]),'Visible','off'); %schemes
set(h2([1]),'Visible','off'); %schemes
COL=colormap(ha(1),'lines');
COL(4,:)=[226 139 138]/255; %afferent color
COL(5,:)=[64 128 0]/255;%[57 181 74]/255; output color
COL(6,:)=[210 148 37]/255;%COL(7,:);%interneuron
COL(7,:)=[86 124 141]/255;%output color

ft=data1.FRAME.t(frameind);
acc(:,1)=data1.FRAME.accx;
acc(:,2)=data1.FRAME.accy;
acc(:,3)=data1.FRAME.accz;
[coef,score,latent,tqs,explained]=pca(acc);
posture=data1.FRAME.posture;
rfin=posture(:,9)/2; rfin=rfin-nanmedian(rfin);
lfin=posture(:,10)/2; lfin=lfin-nanmedian(lfin);
tpos=mean(posture(:,7:8),2);
x1=posture(:,1); x1=(x1-nanmin(x1))/(nanmax(x1)-nanmin(x1));
y1=posture(:,2); y1=(y1-nanmin(y1))/(nanmax(y1)-nanmin(y1));
az=posture(:,3)/pi;
eind=[eod1(seg).ind0(fnum):eod1(seg).ind0(fnum+1)-1];   %eod indices for file
edata=eod1(seg).data(eind,:);   %eod data for file
et=eod1(seg).t(eind,:);    %eod1 time

%% plot ephys data
% ind=find(inrange(t-ft(1),xlim)); %eod indices in t
xlim_=ft(1)+[13.7 14];
ind=find(inrange(t,xlim_));
a1=a(ind);
axes(h2(2));
h2(2).clo;
% h1(4).Position=p1(4,:)+[0 0.025 0 -0.025];
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
Rt=double(data1.RAST(:,1))/ops1.samplerate-data1.FILE.offset;
ind73=find(inrange(Rt,ft([1 end])) & data1.RAST(:,2)==73);
ind=find(ismember(t,Rt(ind73)))'; %eod indices in t
Idx=ind*ones(size(Ix))+ones(size(ind))*Ix;    
Tr73=a(Idx);
H=plot(t(Idx')-xlim_(1),Tr73','Color',COL(6,:),'LineWidth',wline);
ind72=find(inrange(Rt,ft([1 end])) & data1.RAST(:,2)==72);
ind=find(ismember(t,Rt(ind72)))'; %eod indices in t
Idx=ind*ones(size(Ix))+ones(size(ind))*Ix; 
Tr72=a(Idx);
H=plot(t(Idx')-xlim_(1),Tr72','Color',COL(5,:),'LineWidth',wline);
set(gca,'Xlim',[0 0.155],'Ylim',[-500 1000],'FontSize',fontsize);
set(gca,'XTick',[0 0.05 0.1 0.15],'XTickLabel',[0 50 100 150],...
    'YTick',[-500 0 500 1e3],'YTickLabel',{'-0.5','0','0.5','1'},'Box','off')
set(h2(2),'TickLength',[0.025 0.025],'TickDir','both');

% insets
A6=axes();
A6.clo;
dx=1;
set(A6,'Position',h2(2).Position.*[1 1 .55 .5] + [h2(2).Position(3)*.25 h2(2).Position(4)*.5 0 0]);
x=[1:size(eod1(seg).avtraces,2)]/ops1.samplerate*1e3+ops1.traceblnk;
ind=find(inrange(x,[1 6]));
Y=eod1(seg).avtraces(:,:,ch(1),1);
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
% h2(1).Position=p2(1,:)+[0.025 0 -0.025 0];
[hb,pb]=tight_subplot(3, 1, [.0 .0],[.0 .00],[.00 .0],[.2 .5 .3],[],h1(2)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
% [hb,pb]=split_axes(h1(4),2,1,0.0,0.005);
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
set(H,'LineWidth',wline,'Color',COL(1,:));
k=k+1;
H=plot(ft-ft(1),score(frameind,1)+k);
set(H,'LineWidth',wline,'Color',COL(1,:));
k=k+1;
H=plot(ft-ft(1),tpos(frameind)+k);
set(H,'LineWidth',wline,'Color',COL(1,:));
xlim=[0 ft(end)-ft(1)];
set(hb(2),'FontSize',fontsize,'YTick',[],'Xlim',xlim,'Ylim',[-.5 4.5],'XTick',[],'Box','off');

% [hc,pc]=split_axes(hb(1),4,1,0.0,0.005);
isi=diff(data1.EOD.t);
et=data1.EOD.t(inrange(data1.EOD.t,ft([1 end])));
ifr=1./isi(inrange(data1.EOD.t,ft([1 end])));
axes(hb(1));
H=plot(et-ft(1),ifr,'.');
set(H,'MarkerSize',msize,'Color',COL(2,:),'MarkerFaceColor',COL(2,:));
set(hb(1),'FontSize',fontsize,'Xlim',xlim,'Box','off','XTick',[],'Ylim',[4 12],'YTick',[5 10],'TickLength',[0.025 0.25],'TickDir','both');
%%
COL(1,:)=[226 139 138]/255; %afferent color
L0=1;
L1=0.25;
col{1}=colormap_graded(COL(1,:),L0,L1);
COL(2,:)=[210 148 37]/255;%COL(7,:);%interneuron
col{2}=colormap_graded(COL(2,:),L0,L1);
COL(3,:)=[64 128 0]/255;%[57 181 74]/255; output color
col{3}=colormap_graded(COL(3,:),.7,.1);

[h23,p23]=tight_subplot(3, 1, [.0 .0],[.0 .00],[.00 .0],[],[],h2(3)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
msize=1;
axes(h23(1));
h23(1).clo;
H=plot(downsample(Te,K),downsample(U2,K),'.');
set(H,'Color',COL(3,:),'MarkerSize',msize);
hold on;
H=Vline(tswitch(2:end-1));
set(H(:),'LineStyle',':');
set(gca,'XTick',[],'FontSize',fontsize,'Xlim',[T(1)-.1 max(T)],'box','off','Color','none');
set(gca,'TickLength',[tlength tlength],'TickDir','both');

axes(h23(2));
h23(2).clo;
H=plot(downsample(Te,K),downsample(U1,K),'.');
set(H,'Color',COL(2,:),'MarkerSize',msize);
hold on;
H=Vline(tswitch(2:end-1));
set(H(:),'LineStyle',':');
set(gca,'XTick',[],'YTick',[0 2 4],'FontSize',fontsize,'Xlim',[T(1)-.1 max(T)],'box','off','Color','none');
set(gca,'TickLength',[tlength tlength],'TickDir','both');

axes(h23(3));
h23(3).clo;
patch([T;flipud(T)],[Tr(:,5);flipud(Tr(:,1))]/1e3,col{1}(10,:),'FaceAlpha',1,'FaceColor',col{1}(8,:),'EdgeAlpha',0);
hold on;
patch([T;flipud(T)],[Tr(:,4);flipud(Tr(:,2))]/1e3,col{1}(20,:),'FaceAlpha',1,'FaceColor',col{1}(15,:),'EdgeAlpha',0);
H=plot(T,Tr(:,3)/1e3);
set(H,'Color',COL(1,:),'LineWidth',wline);
set(gca,'TickLength',[tlength tlength],'TickDir','both');

H=Vline(tswitch(2:end-1));
set(H(:),'LineStyle',':');
set(gca,'XTick',[0:10:40],'FontSize',fontsize,'Ylim',[.0 .8],'YTick',[.3 .6],'YTickLabelMode','auto','Xlim',[T(1)-.1 max(T)],'Box','off','XTickLabelMode','auto','Color','none');

%% third row - snapshots
% spike data
segs=[3 4 4];
ind_files=[30 2 38];
% get data
chns=find(ops2.outchans);
ch=chns([5,6,7,8]);
% ch2=chns([7,8]);
ch1=[1 2];
ch2=[3 4];
Ix=[-15:1:45];
tsp=Ix/30;
for s=1:numel(segs)
    fname=ops2.seg(segs(s)).files(ind_files(s));
    sesspath=[ops2.datapath,'\',num2str(ops2.seg(segs(s)).dates),'\'];
    [t_ex,amp_ex]=read_bonsai_binary([sesspath,'amp_',num2str(fname)],ops2.samplerate,ops2.chan_num,ops2.blocksize,ch,'adc');
    t_ex=t_ex-file2(segs(s)).offset(ind_files(s));
    for i=1:size(amp_ex,2)
        aa{i}=amp_ex(:,i)-medfilt1(amp_ex(:,i),60);
    end    
    datafile=[sesspath,'data_',num2str(fname)];
    load(datafile);      
    Rt=double(data.RAST(:,1))/ops2.samplerate-data.FILE.offset;
    if(unit1_name(s))
        ind_unit1=find(data.RAST(:,2)==unit1_name(s));    
        ind_unit1=ind_unit1(inrange(ind_unit1,[15 size(amp_ex,1)-45]));
        ind=find(ismember(t_ex,Rt(ind_unit1)))'; %eod indices in t_ex        
        Idx=ind*ones(size(Ix))+ones(size(ind))*Ix;    
        for j=1:numel(ch1)
            Tr1=aa{ch1(j)}(Idx);
            Tr1_q{s,j}=quantile(Tr1,Q)';
        end
    end
    if(unit2_name(s))
        ind_unit2=find(data.RAST(:,2)==unit2_name(s));    
        ind_unit2=ind_unit2(inrange(ind_unit2,[15 size(amp_ex,1)-45]));
        ind=find(ismember(t_ex,Rt(ind_unit2)))'; %eod indices in t        
        Idx=ind*ones(size(Ix))+ones(size(ind))*Ix;    
        for j=1:numel(ch2)
            Tr2=aa{ch2(j)}(Idx);
            Tr2_q{s,j}=quantile(Tr2,Q)';
        end
    end    
end

[h3,p3]=tight_subplot(3, 3, [.01 .05],[.0 .00],[.00 .0],[],[.25 .375 .375],ha(3)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]

% lfp examples
x=ops2.traceblnk+[1:size(eod2(3).avtraces,2)]/30;
indx=find(inrange(x,[0 7]));

for s=1:numel(segs)
    [h31(s,:),p2_(s,:)]=tight_subplot(1,2,[0 0],[0 0],[0 0],[],[],h3(3*s-2));
    for j=1:2
        axes(h31(s,j));
        h31(s,j).clo;
        Atr=sort(eod2(segs(s)).avtraces(:,indx,lfpind(j),ind_files(s)))/1e3;
        x=x(indx);
        P=patch([x fliplr(x)],[Atr(10,:) fliplr(Atr(1,:))],COL(1,:),'FaceAlpha',1,'FaceColor',col{1}(8,:),'EdgeAlpha',0);
        hold on;
        P=patch([x fliplr(x)],[Atr(7,:) fliplr(Atr(3,:))],COL(1,:),'FaceAlpha',1,'FaceColor',col{1}(15,:),'EdgeAlpha',0);    
        H=plot(x,Atr(5,:)');
        set(H,'Color',COL(1,:),'LineWidth',wline);

        if(j==1)
            set(gca,'Ylim',[-.3 .3],'XAxisLocation','bottom','XTick',[2 4 6],'XTickLabel',[],'YTick',[-.2 0 .2],'YTickLabelMode','auto','FontSize',fontsize,'Box','off');
        else
            set(gca,'Ylim',[-.3 .3],'XAxisLocation','bottom','XTick',[2 4 6],'XTickLabel',[],'YTick',[],'YColor','none','FontSize',fontsize,'Box','off');
        end
        set(gca,'TickLength',[tlength tlength],'TickDir','both');
    end
end
set(h31(end,:),'XTick',[2 4 6],'XTickLabelMode','auto');

% unit1 examples
lines=1e2;
binsize=.1;
for s=1:numel(segs)
    [h32(s,:),p3_(s,:)]=tight_subplot(1,3,[0 0],[0 0],[0 0],[],[.2 .2 .6],h3(3*s-1));
    for j=1:2
        axes(h32(s,j));
        h32(s,j).clo;
        P=patch([tsp fliplr(tsp)],[Tr1_q{s,j}(:,5)' fliplr(Tr1_q{s,j}(:,1)')]/1e3,COL(2,:),'FaceAlpha',1,'FaceColor',col{2}(10,:),'EdgeAlpha',0);
        hold on;
        P=patch([tsp fliplr(tsp)],[Tr1_q{s,j}(:,4)' fliplr(Tr1_q{s,j}(:,2)')]/1e3,COL(2,:),'FaceAlpha',1,'FaceColor',col{2}(15,:),'EdgeAlpha',0);        
        H=plot(tsp,Tr1_q{s,j}(:,3)/1e3');        
        set(H,'Color',COL(2,:),'LineWidth',wline);        
        if(j==1)
            set(gca,'Ylim',[-.2 .1],'XAxisLocation','bottom','XTick',[0 1],'XTickLabel',[],'YTick',[-.1 0],'YTickLabelMode','auto','FontSize',fontsize,'Box','off');
        else
            set(gca,'Ylim',[-.2 .1],'XAxisLocation','bottom','XTick',[],'YTick',[],'YColor','none','FontSize',fontsize,'Box','off');
        end
        set(gca,'TickLength',[tlength tlength],'TickDir','both');

    end

    %raster,psth
    axes(h32(s,3));
    h32(s,3).clo;
    h32(s,3).Position=p3_{s,3}+[.05 0 -0.05 0];
    ind=[eod2(segs(s)).ind0(ind_files(s)):eod2(segs(s)).ind0(ind_files(s)+1)-1];
    [ha,hb]=plot_sorted_raster(eod2(segs(s)).raster{unit1col(s)-18},eod2(segs(s)).data(:,18),ops2,ind,lines,1,binsize,col(2),h32(s,3));
    delete(ha(2));
    hb(1).Position(3)=h32(s,3).Position(3);
    hb(2).Position(3)=h32(s,3).Position(3);    
    set(hb(2),'YLim',[0 300],'YTick',[0 200],'XTick',[]);
    
end
set(hb(2),'YLim',[0 300],'YTick',[0 200],'XTick',[-20 0 20 40]);
set(h32(end,:),'XTick',[0 1],'XTickLabelMode','auto');

% unit2 examples
for s=1:numel(segs)
    if(~unit2col(s))
        h3(3*s).Visible='off';
        continue;
    end
    [h33(s,:),p4_(s,:)]=tight_subplot(1,3,[0 0],[0 0],[0 0],[],[.2 .2 .6],h3(3*s));    
    for j=1:2
        axes(h33(s,j));
        h33(s,j).clo;
        P=patch([tsp fliplr(tsp)],[Tr2_q{s,j}(:,5)' fliplr(Tr2_q{s,j}(:,1)')]/1e3,COL(3,:),'FaceAlpha',1,'FaceColor',col{3}(1,:),'EdgeAlpha',0);
        hold on;
        P=patch([tsp fliplr(tsp)],[Tr2_q{s,j}(:,4)' fliplr(Tr2_q{s,j}(:,2)')]/1e3,COL(3,:),'FaceAlpha',1,'FaceColor',col{3}(25,:),'EdgeAlpha',0);
        H=plot(tsp,Tr2_q{s,j}(:,3)/1e3');
        set(H,'Color',COL(3,:),'LineWidth',wline);
        
        if(j==1)
            set(gca,'Ylim',[-.2 .1],'XAxisLocation','bottom','XTick',[0 1],'XTickLabel',[],'YTick',[-.1 0 ],'YTickLabelMode','auto','FontSize',fontsize,'Box','off');
        else
            set(gca,'Ylim',[-.2 .1],'XAxisLocation','bottom','XTick',[0 1],'XTickLabel',[],'YTick',[],'YColor','none','FontSize',fontsize,'Box','off');
        end
    end
    
    %raster,psth
    axes(h33(s,3));
    h33(s,3).clo;
    h33(s,3).Position=p4_{s,3}+[.05 0 -0.05 0];    
    ind=[eod2(segs(s)).ind0(ind_files(s)):eod2(segs(s)).ind0(ind_files(s)+1)-1];
    [ha,hb]=plot_sorted_raster(eod2(segs(s)).raster{unit2col(s)-18},eod2(segs(s)).data(:,17),ops2,ind,lines,1,binsize,col(3),h33(s,3));
    delete(ha(2));
    hb(1).Position(3)=h33(s,3).Position(3);
    hb(2).Position(3)=h33(s,3).Position(3);
    set(hb(2),'YLim',[0 60],'YTick',[0 30],'XTick',[]);    
end
set(hb(2),'YLim',[0 60],'YTick',[0 30],'XTick',[-20 0 20 40]);
set(h33(end,:),'XTick',[0 1],'XTickLabelMode','auto');

print('-clipboard','-dmeta');
