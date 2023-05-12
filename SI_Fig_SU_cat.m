load('z:\analysis_data\ELL_Database.mat') 
if(~exist('eod2'))
    load([rootdir,'analysis_data',sep,'20190624',sep,'output.mat']);
    eod2=eod;
    file2=file;
    ops2=ops;
end

COL=[236 157 118;210 149 0;74 133 34]/255;

burst_th=0.15;
baseline_th1=7;
baseline_th2=3;
ind_type1=find(units.burst>=burst_th & units.baseline<baseline_th1);
ind_type2=find(units.burst<burst_th & units.baseline>baseline_th2);
bgc=[1 1 1];
fontsize=7;
msize=1;
wline=1;
nline=0.25;
d=0.1;

Fall=figure;
[ha, pa] = tight_subplot(2, 1, [d .00],[d .0],[d .0],[.33 .67]); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in
%% first row - examples
[hb,pb] = tight_subplot(1, 3, [.00 d],[.0 .0],[.0 .0],[1],1/3*[1 1 1],ha(1)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in

%rasters + psth
[hb1,pb1] = tight_subplot(3, 1, [.00 0],[.0 .0],[.0 .0],[.25 .25 .5],[1],hb(1)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in
I=[1 3];

bpre=ops2.eodblankpre;
bpost=ops2.eodblankpost;

%raster
N=1e3;
N0=size(eod2(4).data,1);

for i=1:2
    axes(hb1(i));
    hb1(i).clo;
    rast=eod2(4).raster{I(i)};
    rast(:,1)=rast(:,1)*1e3;
    r=randi(N0,N,1);
    ind=find(ismember(rast(:,2),r));
    q=zeros(N0,1);
    q(r)=[1:N];
    drast=rast(ind,:); 
    drast(:,2)=q(drast(:,2));
    rastplot(drast,'|',COL(i+1,:),msize);  
    hold on;
    A=area([-bpre bpost]*1e3,[N N],'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',0.2);
    set(hb1(i),'Color',bgc,'XColor',1-bgc,'YColor',1-bgc,'XTick',[],'YTick',[],'box','off');
    set(hb1(i),'Xlim',[-10 25],'Ylim',[0 N]);
end

%psth
ksmooth=3;
binsize=.5;

edges=[-ops2.rasterpre*1e3:binsize:(-bpre*1e3)  (bpost*1e3):binsize:ops2.rasterpost*1e3];
bins=edge2bin(edges);        
indblank=find(inrange(bins,[-bpre bpost]*1e3));
axes(hb1(3));
hb1(3).clo;
I=[1 3];
for i=1:2
rast=eod2(4).raster{I(i)};
rast(:,1)=rast(:,1)*1e3;
h=histcounts(rast(:,1),edges)/size(eod2(4).data,1)/binsize*1e3;

hh=smooth(h,ksmooth);
hh(indblank)=nan;
H=plot(bins,hh);
set(H,'Color',COL(i+1,:),'LineWidth',wline);
hold on;
end
A=area([-bpre bpost]*1e3,[1e3 1e3],'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',0.2);
set(hb1(3),'Color',bgc,'XColor',1-bgc,'YColor',1-bgc,'FontSize',fontsize,'box','off','YAxisLocation','left');
set(gca,'Xlim',[-10 25],'Ylim',[0 90],'XTick',[-5 0 10 20],'YTick',[40 80]);

% ISI 
axes(hb(2));
hb(2).clo;
isi_x=file2(4).units.isi.bins;
for i=1:2
    isi_y=file2(4).units.isi.val{1}(:,I(i))*100;
    H=plot(isi_x,isi_y);
    set(H,'LineWidth',wline,'Color',COL(i+1,:));
    hold on;
end
set(gca,'Xlim',[0 50],'Ylim',[0 6],'FontSize',fontsize);
 
% burst
fdate=num2str(ops2.seg(4).dates(1));
dname=[ops2.datapath,'\',fdate,'\'];
fnum=1;
fname=ops2.seg(4).files(fnum); %file to take
datafile=[dname,'data_',num2str(fname)];
if(~exist('data'))
    load(datafile);        %for traces
end
i=1;
uname=[32 34];
axes(hb(3));
hb(3).clo;
H=plot([7 7],[0 1],':k');
set(H,'LineWidth',wline);
hold on;
for i=1:2
    rast=double(data.RAST(data.RAST(:,2)==uname(i),1))/ops2.samplerate;
    isi=diff(rast)*1e3;
    bisi=min(isi(1:end-1),isi(2:end));
    axes(hb(3));
    H=cdfplot(bisi);
    set(H,'LineWidth',wline,'Color',COL(i+1,:));
    hold on;
    %insets
    hi(i)=axes();   
    hi(i).clo;
    axes(hi(i));    
    H=scatter(isi(1:N-1),isi(2:N),msize,COL(i+1,:),'filled');
    set(H,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.1);    
    hold on;
    P1=patch([1 1 7 7 1e3 1e3 1],[1 1e3 1e3 7 7 1 1],[0 0 0]);
    set(P1,'FaceAlpha',.1,'FaceColor',[0 0 0],'LineStyle','none')
    set(hi(i),'XScale','log','YScale','log','Xlim',[1 1e3],'Ylim',[1 1e3],'box','on','XTick',[],'YTick',[]);

end
axes(hb(3));
set(gca,'XScale','log','XGrid','off','YGrid','off','XLabel',[],'YLabel',[],'FontSize',fontsize,'Xlim',[1 1e3]);
set(gca,'XTick',[1 7  100 700],'XTickLabel',{'1','7','10','70','100','700'},'Box','off');

hi(1).Position=pb{3}.*[1 1 0.3 0.35]+[pb{3}(3)*.2 pb{3}(4)*.6 0 0];
uistack(hi(1),'top')
hi(2).Position=pb{3}.*[1 1 0.3 0.35]+[pb{3}(3)*.6 pb{3}(4)*.2 0 0];
uistack(hi(2),'top');
%% second row - stats distribution
[hc,pc] = tight_subplot(1, 2, [.00 d],[.0 .0],[.0 .0],[1],[.5 .5],ha(2)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in
[hc1,pc1] = tight_subplot(2, 2, [.00 .0],[.0 .0],[.0 .0],[.2 .8],[.8 .2],hc(1));
[hc2,pc2] = tight_subplot(2, 2, [.00 .0],[.0 .0],[.0 .0],[.2 .8],[.8 .2],hc(2));
hc1(2).Visible='off';
hc2(2).Visible='off';
%% baseline vs. bursts
xlim=[0 25];
ylim=[0 .9];
axes(hc1(3));
hc1(3).clo;
S=scatter(units.baseline(ind_type1),units.burst(ind_type1),msize*30,COL(2,:),'filled');
set(S,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
hold on;
S=scatter(units.baseline(ind_type2),units.burst(ind_type2),msize*30,COL(3,:),'filled','diamond');
set(S,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
H=plot([0 20],[0 0.6],':k');
set(H,'LineWidth',wline,'Color',0.5*[1 1 1],'LineStyle','-.');
set(hc1(3),'XTick',[0 10 20],'Box','off','FontSize',fontsize,'YTick',[0 0.2 .4 .6 .8],'YLim',ylim,'XLim',xlim);

axes(hc1(1));
hc1(1).clo;
bl_x=linspace(0,25,1e2); bl_bw=0.75;
bl_h1=ksdensity(units.baseline(ind_type1),bl_x,'Bandwidth',bl_bw);
bl_h2=ksdensity(units.baseline(ind_type2),bl_x,'Bandwidth',bl_bw);
area(bl_x,bl_h1,'FaceColor',COL(2,:),'FaceAlpha',0.5,'LineStyle','none');
hold on;
area(bl_x,bl_h2,'FaceColor',COL(3,:),'FaceAlpha',0.5,'LineStyle','none');
set(hc1(1),'XTick',[],'Box','off','FontSize',16,'YTick',[],'XLim',xlim,'Ylim',[0 0.35]);

axes(hc1(4));
hc1(4).clo;
br_x=linspace(0,1,1e2); br_bw=0.05;
br_h1=[ksdensity(units.burst(ind_type1),br_x,'Bandwidth',br_bw) 0 0];
br_h2=[ksdensity(units.burst(ind_type2),br_x,'Bandwidth',br_bw) 0 0];
br_x=[br_x 1 0];
patch(br_h1,br_x,COL(2,:),'FaceAlpha',0.5,'LineStyle','none');
hold on;
patch(br_h2,br_x,COL(3,:),'FaceAlpha',0.5,'LineStyle','none');
set(hc1(4),'XTick',[],'Box','off','FontSize',16,'YTick',[],'YLim',ylim);

%%
xlim=[0 350];
ylim=[0 32];
axes(hc2(3));
hc2(3).clo;
S=scatter(units.max(ind_type1),units.spmode(ind_type1),msize*30,COL(2,:),'filled');
set(S,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
hold on;
S=scatter(units.max(ind_type2),units.spmode(ind_type2),msize*30,COL(3,:),'filled','diamond');
set(S,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
set(hc2(3),'XTick',[0 100 200 300],'Box','off','FontSize',fontsize,'YTick',[0 10 20 30],'YLim',ylim,'XLim',xlim);

axes(hc2(1));
hc2(1).clo;
m_x=linspace(0,350,1e2); m_bw=20;
m_h1=ksdensity(units.max(ind_type1),m_x,'Bandwidth',m_bw);
m_h2=ksdensity(units.max(ind_type2),m_x,'Bandwidth',m_bw);
area(m_x,m_h1,'FaceColor',COL(2,:),'FaceAlpha',0.5,'LineStyle','none');
hold on;
area(m_x,m_h2,'FaceColor',COL(3,:),'FaceAlpha',0.5,'LineStyle','none');
set(hc2(1),'XTick',[],'Box','off','FontSize',16,'YTick',[],'XLim',xlim,'Ylim',[0 0.03]);

axes(hc2(4));
hc2(4).clo;
sp_x=linspace(0,30,1e2); sp_bw=1;
sp_h1=[ksdensity(units.spmode(ind_type1),sp_x,'Bandwidth',sp_bw) 0 0];
sp_h2=[ksdensity(units.spmode(ind_type2),sp_x,'Bandwidth',sp_bw) 0 0];
sp_x=[sp_x 1 0];
patch(sp_h1,sp_x,COL(2,:),'FaceAlpha',0.5,'LineStyle','none');
hold on;
patch(sp_h2,sp_x,COL(3,:),'FaceAlpha',0.5,'LineStyle','none');
set(hc2(4),'XTick',[],'Box','off','FontSize',16,'YTick',[],'YLim',ylim);

set(Fall,'Units','centimeters');
set(Fall,'Position', [2  2  12  9],'Color',[1 1 1]);
