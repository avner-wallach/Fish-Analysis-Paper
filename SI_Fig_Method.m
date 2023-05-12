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
    load([rootdir,'analysis_data',sep,'20200113',sep,'output.mat']);
    eod1=eod;
    file1=file;
    ops1=ops;
%     frame1=frame;
end

load('Z:\mormyrid_data\fish_images\fish_silhouette.mat');
str_image=IMGS;
setenv('DATAPATH','Z:\mormyrid_data');
%% plot LFP distribution vs experiment time
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
    for j=1:numel(eod1(segs(i)).ind0)-1
        ind=[eod1(segs(i)).ind0(j):eod1(segs(i)).ind0(j+1)-1];
        lfp=eod1(segs(i)).data(ind,lfpcol)*eod1(segs(i)).LFPs(scaleind(i))+eod1(segs(i)).LFPm(scaleind(i));
%         lfp_q(j,:)=nanmean(lfp)+nanstd(lfp)*[-1 0 1];
        lfp_q(j,:)=quantile(lfp,Q);
%         unit1_q(j,:)=nanmean(unit1)+nanstd(unit1)*[-1 0 1];;        
    end
    isi=eod1(segs(i)).data(:,1);
    isi(isnan(isi))=nanmean(isi);        
    if(unit1col(i))
        unit1=smooth(eod1(segs(i)).data(:,unit1col(i)),K);
    else
        unit1=zeros(size(isi));
    end
    if(unit2col(i))
        unit2=smooth(eod1(segs(i)).data(:,unit2col(i)),K);
    else
        unit2=zeros(size(isi));
    end
%     Tr=cat(3,Tr,permute(eod1(i).avtraces(:,:,lfpind,:),[1 2 4 3]));
    Tr=[Tr;lfp_q];
    U1=[U1;unit1];
    U2=[U2;unit2];
    ISI=[ISI;isi];
    dt=[dt;eod1(segs(i)).t(eod1(segs(i)).ind0(2:end)-1)/3600];
    dtswitch=[dtswitch;sum(eod1(segs(i)).t(eod1(segs(i)).ind0(2:end)-1)/3600)];
end
T=cumsum(dt);
Te=cumsum(ISI)/3600;
tswitch=[1;cumsum(dtswitch)];
%% spike data
segs=[3 4 4];
ind_files=[30 2 38];
% get data
chns=find(ops1.outchans);
ch=chns([5,6,7,8]);
% ch2=chns([7,8]);
ch1=[1 2];
ch2=[3 4];
Ix=[-15:1:45];
tsp=Ix/30;
for s=1:numel(segs)
    fname=ops1.seg(segs(s)).files(ind_files(s));
    sesspath=[ops1.datapath,'\',num2str(ops1.seg(segs(s)).dates),'\'];
    [t,amp]=read_bonsai_binary([sesspath,'amp_',num2str(fname)],ops1.samplerate,ops1.chan_num,ops1.blocksize,ch,'adc');
    t=t-file1(segs(s)).offset(ind_files(s));
    for i=1:size(amp,2)
        a{i}=amp(:,i)-medfilt1(amp(:,i),60);
    end    
    datafile=[sesspath,'data_',num2str(fname)];
    load(datafile);      
    Rt=double(data.RAST(:,1))/ops1.samplerate-data.FILE.offset;
    if(unit1_name(s))
        ind_unit1=find(data.RAST(:,2)==unit1_name(s));    
        ind_unit1=ind_unit1(inrange(ind_unit1,[15 size(amp,1)-45]));
        ind=find(ismember(t,Rt(ind_unit1)))'; %eod indices in t        
        Idx=ind*ones(size(Ix))+ones(size(ind))*Ix;    
        for j=1:numel(ch1)
            Tr1=a{ch1(j)}(Idx);
            Tr1_q{s,j}=quantile(Tr1,Q)';
        end
    end
    if(unit2_name(s))
        ind_unit2=find(data.RAST(:,2)==unit2_name(s));    
        ind_unit2=ind_unit2(inrange(ind_unit2,[15 size(amp,1)-45]));
        ind=find(ismember(t,Rt(ind_unit2)))'; %eod indices in t        
        Idx=ind*ones(size(Ix))+ones(size(ind))*Ix;    
        for j=1:numel(ch2)
            Tr2=a{ch2(j)}(Idx);
            Tr2_q{s,j}=quantile(Tr2,Q)';
        end
    end    
end

%% figure
Fall=figure;
set(Fall,'Units','centimeters');
set(Fall,'Position', [2  2  12  10],'Color',[1 1 1]);
bgc=[0 0 0];
% bgc=[1 1 1];
fontsize=7;
nline=0.25;
wline=.5;
msize=3;
tlength=0.0125;
COL(1,:)=[226 139 138]/255; %afferent color
L0=1;
L1=0.25;
col{1}=colormap_graded(COL(1,:),L0,L1);
COL(2,:)=[210 148 37]/255;%COL(7,:);%interneuron
col{2}=colormap_graded(COL(2,:),L0,L1);
COL(3,:)=[64 128 0]/255;%[57 181 74]/255; output color
col{3}=colormap_graded(COL(3,:),.7,.1);

[ha, pos] = tight_subplot(2, 1, [.1 .0],[.05 .05],[.1 .025],[]); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
%% first row: implant, long term dynamics
ha(1).Position=pos{1}+[0.2 0.0 -0.2 0];
[h1,p1]=tight_subplot(3, 1, [.0 .0],[.0 .00],[.00 .0],[],[],ha(1)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]

axes(h1(1));
h1(1).clo;
H=plot(downsample(Te,K),downsample(U2,K),'.');
set(H,'Color',COL(3,:),'MarkerSize',msize);
hold on;
H=Vline(tswitch(2:end-1));
set(H(:),'LineStyle',':');
set(gca,'XTick',[],'FontSize',fontsize,'Xlim',[T(1)-.1 max(T)],'box','off');
set(gca,'TickLength',[tlength tlength],'TickDir','both');

axes(h1(2));
h1(2).clo;
H=plot(downsample(Te,K),downsample(U1,K),'.');
set(H,'Color',COL(2,:),'MarkerSize',msize);
hold on;
H=Vline(tswitch(2:end-1));
set(H(:),'LineStyle',':');
set(gca,'XTick',[],'YTick',[0 2 4],'FontSize',fontsize,'Xlim',[T(1)-.1 max(T)],'box','off');
set(gca,'TickLength',[tlength tlength],'TickDir','both');

axes(h1(3));
h1(3).clo;
patch([T;flipud(T)],[Tr(:,5);flipud(Tr(:,1))]/1e3,col{1}(10,:),'FaceAlpha',1,'FaceColor',col{1}(8,:),'EdgeAlpha',0);
hold on;
patch([T;flipud(T)],[Tr(:,4);flipud(Tr(:,2))]/1e3,col{1}(20,:),'FaceAlpha',1,'FaceColor',col{1}(15,:),'EdgeAlpha',0);
H=plot(T,Tr(:,3)/1e3);
set(H,'Color',COL(1,:),'LineWidth',wline);
set(gca,'TickLength',[tlength tlength],'TickDir','both');

H=Vline(tswitch(2:end-1));
set(H(:),'LineStyle',':');
set(gca,'XTick',[0:10:40],'FontSize',fontsize,'Ylim',[.2 .7],'YTick',[.2 .4 .6],'YTickLabelMode','auto','Xlim',[T(1)-.1 max(T)],'Box','off','XTickLabelMode','auto');

%% swcond row - snapshots
[h2,p2]=tight_subplot(3, 3, [.1 .05],[.0 .00],[.00 .0],[],[.25 .375 .375],ha(2)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]

% lfp examples
x=ops1.traceblnk+[1:size(eod1(3).avtraces,2)]/30;
indx=find(inrange(x,[0 7]));

for s=1:numel(segs)
    [h2_(s,:),p2_(s,:)]=tight_subplot(1,2,[0 0],[0 0],[0 0],[],[],h2(3*s-2));
    for j=1:2
        axes(h2_(s,j));
        h2_(s,j).clo;
        Atr=sort(eod1(segs(s)).avtraces(:,indx,lfpind(j),ind_files(s)))/1e3;
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
set(h2_(end,:),'XTick',[2 4 6],'XTickLabelMode','auto');

% unit1 examples
lines=1e2;
binsize=.1;
for s=1:numel(segs)
    [h3_(s,:),p3_(s,:)]=tight_subplot(1,3,[0 0],[0 0],[0 0],[],[.2 .2 .6],h2(3*s-1));
    for j=1:2
        axes(h3_(s,j));
        h3_(s,j).clo;
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
    axes(h3_(s,3));
    h3_(s,3).clo;
    h3_(s,3).Position=p3_{s,3}+[.05 0 -0.05 0];
    ind=[eod1(segs(s)).ind0(ind_files(s)):eod1(segs(s)).ind0(ind_files(s)+1)-1];
    [ha,hb]=plot_sorted_raster(eod1(segs(s)).raster{unit1col(s)-18},eod1(segs(s)).data(:,18),ops1,ind,lines,1,binsize,col(2),h3_(s,3));
    delete(ha(2));
    hb(1).Position(3)=h3_(s,3).Position(3);
    hb(2).Position(3)=h3_(s,3).Position(3);    
    set(hb(2),'YLim',[0 300],'YTick',[0 200]);
    
end
set(h3_(end,:),'XTick',[0 1],'XTickLabelMode','auto');

% unit2 examples
for s=1:numel(segs)
    if(~unit2col(s))
        h2(3*s).Visible='off';
        continue;
    end
    [h4_(s,:),p4_(s,:)]=tight_subplot(1,3,[0 0],[0 0],[0 0],[],[.2 .2 .6],h2(3*s));    
    for j=1:2
        axes(h4_(s,j));
        h4_(s,j).clo;
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
    axes(h4_(s,3));
    h4_(s,3).clo;
    h4_(s,3).Position=p4_{s,3}+[.05 0 -0.05 0];    
    ind=[eod1(segs(s)).ind0(ind_files(s)):eod1(segs(s)).ind0(ind_files(s)+1)-1];
    [ha,hb]=plot_sorted_raster(eod1(segs(s)).raster{unit2col(s)-18},eod1(segs(s)).data(:,17),ops1,ind,lines,1,binsize,col(3),h4_(s,3));
    delete(ha(2));
    hb(1).Position(3)=h4_(s,3).Position(3);
    hb(2).Position(3)=h4_(s,3).Position(3);
    set(hb(2),'YLim',[0 60],'YTick',[0 30]);
    
end
set(h4_(end,:),'XTick',[0 1],'XTickLabelMode','auto');