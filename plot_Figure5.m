% make panels for Figure 1: Afferent input in freely moving fish
[ret,sname]=system('hostname');
if(strfind(sname,'Mac'))
    rootdir='/Volumes/aw3057/';
    sep='/';
else
    rootdir='Z:\';
    sep='\';
end
%% load data
if(~exist('eod1'))
    load([rootdir,'analysis_data',sep,'20190131',sep,'output.mat']);
    eod1=eod;
    file1=file;
    ops1=ops;
%     frame1=frame;
end
if(~exist('eod2'))
    load([rootdir,'analysis_data',sep,'20190624',sep,'output.mat']);
    eod2=eod;
    file2=file;
    ops2=ops;
%     frame2=frame;
end

load('Z:\mormyrid_data\fish_images\fish_silhouette.mat');
str_image=IMGS;
setenv('PXLSIZE','0.423'); %pxl to mm conversion
pxl2mm=0.423;
setenv('FRAMESCALE',num2str(ops.framescale));
setenv('SESSDATE','20190131');
setenv('DATAPATH','Z:\mormyrid_data');
load('ROIs.mat','p_IO');
%% plot ex afference, re afference maps
Cslope=[-1 1];
bgcol=0*[1 1 1];
upsamp=1.5;
minsamp=10;
nline=0.25;
wline=0.5;
fontsize=7;
%%
lfp_eod=eod1(4); lfp_file=file1(4); lfp_col=17; lfp_ops=ops1;%lfp_zscale=[430 40];
lfp_zscale=[.730 .230]; lfp_clim=[.5 1];
% lfp_eod=eod2(6); lfp_file=file2(6); lfp_col=21;
unit1_eod=eod2(6); unit1_file=file2(6); unit1_col=24; unit1_ind=ops2.seg(6).ind_lim(2,:); unit1_clim=[0 80]; unit1_lfpcol=21; unit1_rastind=2;
% unit1_eod=eod1(4); unit1_file=file1(4); unit1_col=19; unit1_ind=ops1.seg(4).ind_lim(1,:); unit1_clim=[0 .65]; unit1_lfpcol=17; unit1_rastind=1;
unit2_eod=eod1(4); unit2_file=file1(4); unit2_col=20; unit2_ind=ops1.seg(4).ind_lim(2,:); unit2_clim=[0 80]; unit2_lfpcol=17; unit2_rastind=2;
% unit2_eod=eod2(6); unit2_file=file2(6); unit2_col=23; unit2_ind=ops2.seg(6).ind_lim(1,:); unit2_clim=[0 1.5]; unit2_lfpcol=21; unit2_rast=unit2_eod.raster{1};
xlim=[-550 550];
ylim=[-550 250];
unit1_cmap=ones(3,1)*COL(2,:);
unit2_cmap=ones(3,1)*COL(3,:);
N=4;
%% rois
rois=0;
if(rois)    
    Froi=figure;
    a1=subplot(1,3,1); a2=subplot(1,3,2); a3=subplot(1,3,3);     
    
    axes(a2);
    a2.clo;
    [M1,N1,s1]=plot_tuning2d(unit1_eod,unit1_file,unit1_col,'mfunc','mean','image',str_image.cdata,...
        'clim',unit1_clim,'bgcol',bgcol,'x_lim',xlim,'y_lim',ylim,'y_nbins',y_nbins,...
        'minsamp',minsamp,'upsamp',upsamp,'zscale',[0 40],'ind_lim',unit1_ind,...
        'rastind',unit1_rastind,'ops',ops1,'cmaps',unit1_cmap,'psth_lines',1);
    
    axes(a3);
    a3.clo;
    [M2,N2,s2]=plot_tuning2d(unit2_eod,unit2_file,unit2_col,'mfunc','mean','image',str_image.cdata,...
        'clim',unit2_clim,'bgcol',bgcol,'x_lim',xlim,'y_lim',[-550 300],'y_nbins',y_nbins,...
        'minsamp',minsamp,'upsamp',upsamp,'zscale',[0 40],'ind_lim',unit2_ind,...
        'rastind',unit2_rastind,'ops',ops2,'cmaps',unit2_cmap,'psth_lines',1,'roinum',1);
    
    axes(a1);
    a1.clo;
    [M0,N0,s0]=plot_tuning2d(lfp_eod,lfp_file,lfp_col,'mfunc','mean','image',str_image.cdata,...
        'clim',Coff,'bgcol',bgcol,'x_lim',xlim,'y_lim',ylim,...
        'minsamp',minsamp,'upsamp',upsamp,'ops',lfp_ops,'roinum',N,'cmaps',lfp_cmap);
    p_brass=s0.p;
    
    axes(a2);
    a2.clo;
    [M1,N1,s1]=plot_tuning2d(unit1_eod,unit1_file,unit1_col,'mfunc','mean','image',str_image.cdata,...
        'clim',unit1_clim,'bgcol',bgcol,'x_lim',xlim,'y_lim',ylim,'y_nbins',y_nbins,...
        'minsamp',minsamp,'upsamp',upsamp,'zscore',0,'ind_lim',unit1_ind,...
        'rastind',unit1_rastind,'ops',ops1,'poly',p_brass,'cmaps',unit1_cmap,'psth_lines',1);
    
    axes(a3);
    a3.clo;
    [M2,N2,s2]=plot_tuning2d(unit2_eod,unit2_file,unit2_col,'mfunc','mean','image',str_image.cdata,...
        'clim',unit2_clim,'bgcol',bgcol,'x_lim',xlim,'y_lim',ylim,'y_nbins',y_nbins,...
        'minsamp',minsamp,'upsamp',upsamp,'zscore',0,'ind_lim',unit2_ind,...
        'rastind',unit2_rastind,'ops',ops2,'poly',p_brass,'cmaps',unit2_cmap,'psth_lines',1);
    
    close(Froi);
    save('scatter_polygons.mat','p_brass','-append');
    
else
    load('scatter_polygons.mat');
end

%%
Fall=figure;
set(Fall,'Units','centimeters');
set(Fall,'OuterPosition', [2  2  11.4  15]);
[h0, p0] = tight_subplot(3, 1, [.035 .0],[.025 0.0],[.0 0.0],[.45 .275 .275],[]); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
%% 1st row
[h1, p1] = tight_subplot(1, 2, [.0 .0],[.0 0.0],[0 0.0],[],[.4 .6],h0(1)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
%% get snapshot images
getimages=0;
if(getimages)
    seg=6; fnum=9;
    fdate=num2str(ops2.seg(seg).dates(1));
    dname=[ops2.datapath,'\',fdate,'\'];
    fname=ops2.seg(seg).files(fnum); %file to take
    datafile=[dname,'data_',num2str(fname)];
    load(datafile);        %for traces
    x=data.FRAME.posture(:,1);
    y=data.FRAME.posture(:,2);
    a=data.FRAME.posture(:,3);
    
    [r,th]=convert_exo_polar(x,y,a,file2(seg).objects,pxl2mm);
    % ind=find(inrange(r,[49 51]) & inrange(th,pi/2+[-pi/20 pi/20]));
    x_=r.*sin(th); y_=r.*cos(th);
%     t=frame2(seg).data(:,9);
    Fre=figure;
    sX=size(str_image.cdata,2);
    sY=size(str_image.cdata,1);
    ix=([1:sX]-sX/2)*pxl2mm;
    iy=(-[1:sY]+sY/3)*pxl2mm;
    imagesc(ix,iy,str_image.cdata);  
    hold on;
    plot(x_,y_,'.');
    axis([-175 175 -200 75]);
    set(gca,'YDir','normal');
    for i=1:3
        K=30;
        figure(Fre);
        h=impoly();
        p=h.getPosition;
        ind = find(inpolygon(x_,y_,p(:,1),p(:,2)));
        %get indices for frames
%         [fileidx{i},frameidx{i}]=get_frame_indices(ind,frame(2).ind0);
%         T=mat2cell(t(ind),cellfun(@numel,frameidx{i}));
%         [m,im]=max(cellfun(@numel,frameidx{i}));
        [as,is]=sort(a(ind));
        frind{i}=ind(is(floor(linspace(1,numel(is),K))));        
        figure;
        [frames{i},ax_grid] = get_fish_image(ops2.seg(seg).files(fnum),frind{i}(1:30),'headfix','off');
        close(gcf);
    end
    save('fish_frames_brass','frames','frind','-v7.3');
else
    if(~exist('frames'));
        load('fish_frames_brass','frames','frind');
    end
end

%% get joint figure of all fish
F(:,:,:,1)=frames{1}(10).cdata;
F(:,:,:,2)=frames{2}(20).cdata;
F(:,:,:,3)=frames{3}(1).cdata;
F(:,:,:,4)=frames{4}(24).cdata;
F(:,:,:,5)=frames{5}(6).cdata;
D=double(F);
th=45;
M=median(D,4);
x=1:size(M,2); y=1:size(M,1);
[X,Y]=meshgrid(x,y);
x0=432; y0=400; r=size(M,1)/2;
B0=(((X-x0).^2+(Y-y0).^2)<r.^2);
col=brewermap(9,'Pastel1');
for i=1:5
    B(:,:,:,i)=(abs(D(:,:,:,i)-M)>th);
    c(:,:,1)=col(i,1)*ones(size(M(:,:,1)));
    c(:,:,2)=col(i,2)*ones(size(M(:,:,1)));
    c(:,:,3)=col(i,3)*ones(size(M(:,:,1)));
    FF(:,:,:,i)=(255-D(:,:,:,i)).*c;
end
% figure;
axes(h1(1));
h1(1).clo;
I=image(uint8(M));
set(I,'AlphaData',B0);
hold on;
for i=1:5
    I=image(uint8(FF(:,:,:,i)));
    set(I,'AlphaData',B(:,:,1,i).*B0);
end
set(gca,'XColor','none','YColor','none');
axis('image');
%% exmaple stats
[h12, p12] = tight_subplot(N, 3, [.02 .01],[.0 0.0],[.15 0.0],[],[],h1(2)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]

% h1(2).Visible='off';
%% 2nd,3rd row
y_nbins=25;
COL=[236 157 118;210 149 0;74 133 34]/255;
lfp_cmap=ones(4,1)*COL(1,:);
unit1_cmap=ones(4,1)*COL(2,:);
unit2_cmap=ones(4,1)*COL(3,:);

[h2, p2] = tight_subplot(1, 3, [.0 .0],[.0 0.0],[.025 0.0],[],[],h0(2)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]

% Coff=[-1.5 1.5];
upsamp=2;
minsamp=20;
mfunc='mean';

axes(h2(1));
h2(1).clo;
[M0,N0,s0]=plot_tuning2d(lfp_eod,lfp_file,lfp_col,'mfunc',mfunc,'image',str_image.cdata,...
    'clim',lfp_clim,'bgcol',bgcol,'x_lim',xlim,'y_lim',ylim,'zscale',lfp_zscale,...
    'minsamp',minsamp,'upsamp',upsamp,'ops',lfp_ops,'poly',p_brass,'cmaps',lfp_cmap);
inds=s0.inds;

Q=linspace(0,1,10);
trnum=lfp_ops.seg(4).LFPgroups{4};%str2num(lfp_eod.fnames{lfp_col}(end))};
avtraces=lfp_eod.avtraces(:,:,trnum,1);
J=2;
X=lfp_ops.traceblnk+[1:420]/30;
x=linspace(-2.5,2.5,1e2);

for i=1:N
%     [h12(i,:), p2_(i,:)] = tight_subplot(3, 2, [.02 .02],[.0 0.0],[.0 0.0],[],[],h2(i)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
    m=nanmean(lfp_eod.data(inds{i},lfp_col));
    q1=quantile(lfp_eod.data(inds{i},lfp_col),.25);
    q3=quantile(lfp_eod.data(inds{i},lfp_col),.75);
%     s=nanstd(lfp_eod.data(inds{i},lfp_col));
    q=comp_percentile(lfp_eod.data(:,lfp_col),[q1 m q3]);
    Y=sort(interp1(Q,avtraces(:,:,J),q));
    
    axes(h12(i*3-2));
    h12(i*3-2).clo;
    H=patch([X fliplr(X)],[Y(1,:) fliplr(Y(3,:))],COL(1,:));
    set(H,'LineStyle','none','FaceAlpha',.5);
    hold on;
    H=plot(X,Y(2,:));
    set(H,'Color',COL(1,:),'LineWidth',wline);
    T=text(6,-200,[num2str(round(m*lfp_zscale(2)*1e3+lfp_zscale(1)*1e3)),'\muV'],'FontSize',5,'FontName','Arial','Color',COL(1,:));
    set(gca,'XColor','none','YColor','none','Ylim',[-1e3 400],'Xlim',[1 8],'Color','none');
    
%     axes(h2_(i,2));    
%     [h,x]=ksdensity(lfp_eod.data(inds{i},lfp_col),x);
%     H=plot(x,h);
%     set(H,'LineWidth',wline,'Color',COL(1,:));
%     set(gca,'XAxisLocation','origin','YAxisLocation','origin',...
%         'XTick',[],'YTick',[],'Color','none','Ylim',[0 1],'XLim',[x(1) x(end)],'box','off');    
end
x0=5; y0=-800;
H=plot(x0+[0 1],y0+[0 0],'k',x0+[0 0],y0+[0 250],'k');

faux=figure; haux=axes;
axes(h2(2));
h2(2).clo;
[M1,N1,s1]=plot_tuning2d(unit1_eod,unit1_file,unit1_col,'mfunc',mfunc,'image',str_image.cdata,...
    'clim',unit1_clim,'bgcol',bgcol,'x_lim',xlim,'y_lim',ylim,'y_nbins',y_nbins,...
    'minsamp',minsamp,'upsamp',upsamp,'zscale',[0 40],'ind_lim',unit1_ind,...
    'rastind',unit1_rastind,'ops',ops1,'poly',p_brass,'aux_axes',haux,'cmaps',unit1_cmap,'psth_lines',1);

for i=1:N
    [h12_1(i,:), p12_1(i,:)] = tight_subplot(2, 1, [.0 .0],[.0 0.0],[.0 0.0],[.35 .65],[],h12(i*3-1)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
    for j=1:2
        h12_1(i,j).clo;
        copyaxes(s1.b_out(i,j),h12_1(i,j),1);
        set(h12_1(i,j),'Colormap',s1.b_out(i,j).Colormap,'Color','none');        
    end   
    hold on;
    r=nanmean(unit1_eod.data(s1.inds{i},unit1_col))*40;
    T=text(30,110,[num2str(round(r)),'Hz'],'FontSize',5,'FontName','Arial','Color',COL(2,:));

    
%     h12_(i,4).clo;
%     copyaxes(s1.a_out(i+3),h2_(i,4),1);
%     set(h2_(i,4),'Colormap',s1.a_out(3+i).Colormap,'Color','none');    
end        

axes(h2(3));
[M2,N2,s2]=plot_tuning2d(unit2_eod,unit2_file,unit2_col,'mfunc',mfunc,'image',str_image.cdata,...
    'clim',unit2_clim,'bgcol',bgcol,'x_lim',xlim,'y_lim',ylim,'y_nbins',y_nbins,...
    'minsamp',minsamp,'upsamp',upsamp,'zscale',[0 40],'ind_lim',unit2_ind,...
    'rastind',unit2_rastind,'ops',ops2,'poly',p_brass,'aux_axes',haux,'cmaps',unit2_cmap,'psth_lines',1);

for i=1:N
    [h12_2(i,:), p12_2(i,:)] = tight_subplot(2, 1, [.0 .0],[.0 0.0],[.0 0.0],[.35 .65],[],h12(i*3)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
    for j=1:2
        h12_2(i,j).clo;
        copyaxes(s2.b_out(i,j),h12_2(i,j),1);
        set(h12_2(i,j),'Colormap',s2.b_out(i,j).Colormap,'Color','none');
    end
    hold on;
    r=nanmean(unit2_eod.data(s2.inds{i},unit2_col))*40;    
    T=text(30,110,[num2str(round(r)),'Hz'],'FontSize',5,'FontName','Arial','Color',COL(3,:));
    
%     h2_(i,6).clo;
%     copyaxes(s2.a_out(i+3),h2_(i,6),1);
%     set(h2_(i,6),'Colormap',s2.a_out(3+i).Colormap,'Color','none');
    
end        
aaa=findobj('Type','Area');
set(aaa(:),'ShowBaseline','off');
set(h12_2(:,2),'Ylim',[0 180]);
set(h12_1(:,2),'Ylim',[0 180]);

for i=1:numel(p_brass)
    PP(i,:)=[min(p_brass{i}) range(p_brass{i})];
    rectangle(h2(1),'Position',PP(i,:),'Curvature',[1 1],'EdgeColor',[0 0 0 0.5]);
%     rectangle(h2(2),'Position',P(i,:),'Curvature',[1 1]);
%     rectangle(h2(3),'Position',P(i,:),'Curvature',[1 1]);
end
C=findobj('Type','Colorbar');
set(C(:),'Location','west','AxisLocation','out','LineWidth',nline,'FontSize',fontsize,'TickLength',.02,'TickDirection','both');
for i=1:numel(C)
%     cp(i,:)=C(i).Position;
    C(i).Position(3)=cp(i,3)-0.01;
    C(i).Position(1)=cp(i,1)-0.01;
end
%% third row
[h3, p3] = tight_subplot(1, 2, [.0 .1],[.0 0.0],[.075 0.075],[],[.6 .4],h0(3)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]

m0=discretize(M0(~isnan(M0)),50);
m1=discretize(M1(~isnan(M1)),50);
m2=discretize(M2(~isnan(M2)),50);

BW=1;
f(1,:)=ksdensity(m0,[1:50],'Bandwidth',BW);
f(2,:)=ksdensity(m1,[1:50],'Bandwidth',BW);
f(3,:)=ksdensity(m2,[1:50],'Bandwidth',BW);
[h(1),s(1)]=entropy(m0);
[h(2),s(2)]=entropy(m1);
[h(3),s(3)]=entropy(m2);

x=linspace(0,1,50);
h31=axes();
set(h31,'Position',p3{1}.*[1 1 .15 .4] + [p3{1}(3)*.8 p3{1}(4)/2 0 0]);

fp=3;
h3(1).clo;
h31.clo;
for i=1:3
    axes(h3(1));
    P=patch([i-f(i,:)*fp i+fliplr(f(i,:)*fp)],[x fliplr(x)],COL(i,:));
    set(P,'FaceColor',COL(i,:),'LineStyle','none');
%     H=plot(x,f(i,:));
%     set(H,'Color',COL(i,:),'LineWidth',wline);
    hold on;
    
    axes(h31);
    yyaxis right;
    H=bar(i,s(i));
    set(H,'FaceColor',COL(i,:),'LineStyle','none');
    hold on;
end
set(h3(1),'FontSize',fontsize,'Color','none','Box','off','Xlim',[.5 4],'XTick',[],'YTickLabelMode','auto');
set(h31,'Color','none','XTick',[],'YTick',[.1 .2 .3],'YLim',[0 s(3)*1.1],'FontSize',5,'box','on','TickDir','both','TickLength',[.05 .05],'YColor',[0 0 0]);
yyaxis left;
set(h31,'YTick',[],'YColor',[0 0 0],'Color','none');

%% stats

load('z:\analysis_data\ELL_Database.mat')  
burst_th=0.15;
baseline_th1=7;
baseline_th2=3;
ind_type1=find(units.burst>=burst_th & units.baseline<baseline_th1);
ind_type2=find(units.burst<burst_th & units.baseline>baseline_th2);
p=1.5;
p0=[0.005 0 -0.055 0];
%spatial selectivity
x1=lfp.Sx;
x2=units.Sx(ind_type1);
x3=units.Sx(ind_type2);
X=[x1; x2; x3];
G=[ones(size(x1)); 2*ones(size(x2)); 3*ones(size(x3))];
axes(h3(2));
h3(2).clo;
boxplot(X,G,'PlotStyle','compact','Colors',COL,'Symbol','.w');
hold on;
H=scatter(.05*randn(size(x1))+1,x1,5,'filled');
set(H,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0,'MarkerFaceColor',lighter(COL(1,:),-.5));
H=scatter(.05*randn(size(x2))+2,x2,5,'filled');
set(H,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0,'MarkerFaceColor',lighter(COL(2,:),-.5));
H=scatter(.05*randn(size(x3))+3,x3,5,'filled');
set(H,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0,'MarkerFaceColor',lighter(COL(3,:),-.5));

set(gca,'YLim',[0.0 max(X)*1.15],'YTickMode','auto','YTickLabelMode','auto','XTickLabel',[],'Xlim',[0.5 3.5]);
set(gca,'TickLength',[.025 .025],'TickDir','both');
set(gca,'Color','none','XColor',1-bgcol,'YColor',1-bgcol,'FontSize',fontsize,'box','off');
Sx_p12=boottest(x1,x2,1e5)
Sx_p13=boottest(x1,x3,1e5)
Sx_p23=boottest(x2,x3,1e5)

bwidth=10;
bx=findobj('Tag','Box');
set(bx,'LineWidth',bwidth);
ws=findobj('Tag','Whisker');
set(ws,'LineWidth',wline);
cro=findobj('Tag','MedianInner');
set(cro,'Marker','none');
cri=findobj('Tag','MedianOuter');
set(cri,'MarkerFaceColor',[1 1 1],'MarkerSize',3);
ol=findobj('Tag','Outliers');
delete(ol); %remove outliers
set(h3(2),'Position',p3{2});
%%
set(Fall,'Color','none');
hall=findobj();
% for i=2:numel(hall)
%     hall(i).Visible='on';
% end

srf=findobj('Type','Surface');
img=findobj('Type','Image');
set(srf,'Visible','off');
set(img,'Visible','off');
% export to metafile
print('-clipboard','-dmeta');

for i=2:numel(hall)
    if(~strcmp(hall(i).Type,'root') & ~strcmp(hall(i).Type,'figure'))
        hall(i).Visible='off';
    end
end
set(srf,'Visible','on');
set(img,'Visible','on');
print('-clipboard','-dbitmap','-r720');

