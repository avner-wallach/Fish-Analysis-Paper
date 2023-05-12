% plot acute data figure

%% import data into single struct
getdata=0;
if(getdata)
    strct.pairA_E=[];
    strct.pairA_I=[];
    strct.pairB_E=[];
    strct.pairB_I=[];

    strct.probeA_E=[];
    strct.probeA_I=[];
    strct.probeB_E=[];
    strct.probeB_I=[];

    strct.probeC_E=[];
    strct.probeC_I=[];

    strct.indA_E=[];
    strct.indA_I=[];
    strct.indC_E=[];
    strct.indC_I=[];
    strct.t=[];

    database.nodelay=strct;
    database.ctl=strct;
    database.chin=strct;
    database.trunk=strct;

    cd('D:\contextpairingdata_nbs\');
    D1=dir;
    for d1=3:numel(D1) %go over all experiments
        if(D1(d1).isdir)
            expname=D1(d1).name;
            D2=dir(expname);
            for d2=3:numel(D2) % go over all segments
                dname=D2(d2).name;
                if(strfind(dname,'ctrl'))
                    db='ctl';
                elseif(strfind(dname,'nodelay'))
                    db='nodelay';
                elseif(strfind(dname,'chin'))
                    db='chin';
                else
                    db='trunk';
                end
                V=[];
                %getdata
                D3=what([expname filesep dname]);
                if(numel(D3.mat)==0)
                    continue;
                end
                for d3=1:numel(D3.mat)
                    S=load([expname filesep dname filesep D3.mat{d3}]);
                    sf=fieldnames(S);
                    t=S.(sf{1}).times(1:150);
                    v=S.(sf{1}).values(1:150);
    %                 T=[T t];
                    V=[V v];
                end

                %add data to relevant matrix
                if(strfind(lower(dname),'paira')) %Pair A
                    if(strfind(lower(dname),'ecell')) %E cells
                        database.(db).pairA_E=[database.(db).pairA_E V];
                    else %I cells
                        database.(db).pairA_I=[database.(db).pairA_I V];
                    end
                end
                if(strfind(lower(dname),'pairb')) %Pair B
                    if(strfind(lower(dname),'ecell')) %E cells
                        database.(db).pairB_E=[database.(db).pairB_E V];
                    else %I cells
                        database.(db).pairB_I=[database.(db).pairB_I V];
                    end
                end
                if(strfind(lower(dname),'probee')) %Probe A
                    if(strfind(lower(dname),'ecell')) %E cells
                        database.(db).probeA_E=[database.(db).probeA_E V];
                    else %I cells
                        database.(db).probeA_I=[database.(db).probeA_I V];
                    end
                end
                if(strfind(lower(dname),'probez')) %Probe B
                    if(strfind(lower(dname),'ecell')) %E cells
                        database.(db).probeB_E=[database.(db).probeB_E V];
                    else %I cells
                        database.(db).probeB_I=[database.(db).probeB_I V];
                    end
                end
                if(strfind(lower(dname),'probek')) %Probe C
                    if(strfind(lower(dname),'ecell')) %E cells
                        database.(db).probeC_E=[database.(db).probeC_E V];
                    else %I cells
                        database.(db).probeC_I=[database.(db).probeC_I V];
                    end
                end
                if(numel(strfind(lower(dname),'eind')) | numel(strfind(lower(dname),'inde'))) %ind. A
                    if(strfind(lower(dname),'ecell')) %E cells
                        database.(db).indA_E=[database.(db).indA_E V];
                    else %I cells
                        database.(db).indA_I=[database.(db).indA_I V];
                    end
                end
                if(numel(strfind(lower(dname),'kind')) | numel(strfind(lower(dname),'indk'))) %ind. C
                    if(strfind(lower(dname),'ecell')) %E cells
                        database.(db).indC_E=[database.(db).indC_E V];
                    else %I cells
                        database.(db).indC_I=[database.(db).indC_I V];
                    end
                end
                database.(db).t=t;            
            end
        end
    end

    save('C:\Users\sawte\Dropbox\efish_paper\matlab\acute_data','database');
else
    load('C:\Users\sawte\Dropbox\efish_paper\matlab\acute_data','database');
end
%% generate raster
genrast=1;
if(genrast)
    %get data
    database.EcontRasters=get_rasters('Z:\contextpairingdata_nbs\20211229_011_rasterecell_0_175.mat',0);
    database.IcontRasters=get_rasters('Z:\contextpairingdata_nbs\20211215_020raster_icell_0_200.mat',0);
%      database.Erasters=get_rasters('D:\contextpairingdata_nbs\20211215_raster_ecell_570_895.mat',0);
%      database.ECrasters=get_rasters('D:\contextpairingdata_nbs\20211215_raster_ecell_probek.mat',1);
%      database.Irasters=get_rasters('D:\contextpairingdata_nbs\20220103_016_raster_icell_0_85.mat',0);
%      database.ICrasters=get_rasters('D:\contextpairingdata_nbs\20220103_016_raster_icell_probek.mat',1);
%      save('C:\Users\sawte\Dropbox\efish_paper\matlab\acute_data','database');
end                        
    
%% figure
% bgc=[0 0 0];
bgc=[1 1 1];
fontsize=7;
nline=0.25;
wline=.5;
msize=3;

Fall=figure;
set(Fall,'Units','centimeters');
set(Fall,'OuterPosition', [2  2  11.4  17],'Color',bgc);
COL=brewermap(10,'Set1');
COL(3,:)=[74 133 34]/255;
COL(4,:)=[112 130 56]/255;
COL1=[13 114 186;226 139 138;128 64 152;236 177 32]/255;
[ha, pos] = tight_subplot(2, 1, [.05 .0],[.05 .05],[.05 .0],[.2 .8],[]); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
dbig=.075;
dsmall=.025;
%% first row - schemes
pwidth=.03;
T=1;
ctimes=[0.2 0.5 0.75];
stimes=ctimes+.01;
rtimes=ctimes+0.01;

t1=linspace(0,T,1e3);
stim=zeros(size(t1));
remote=stim;

for i=1:numel(ctimes)
    stim(inrange(t1-stimes(i),[0 pwidth]))=.5;
    remote(inrange(t1-rtimes(i),[0 pwidth]))=.5;
end

[h1,pos1] = tight_subplot(1,3,[0 2*dsmall],[0 0],[0 0.0],[],[.5 .2 .3],ha(1));
h1(1).Visible='off';
%pairing
[h11,pos11] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h1(2));
axes(h11(1));
h11(1).clo;
H=plot(t1,zeros(size(t1)));
set(H,'Color',COL1(1,:),'LineWidth',wline);
hold on;
H=Vline(0,4,ctimes');
set(H,'Color',COL1(1,:),'LineWidth',2*wline,'LineStyle',':');
H=plot(t1,stim+1);
set(H,'Color',COL1(2,:),'LineWidth',wline);
H=plot(t1,remote+2);
set(H,'Color',COL1(3,:),'LineWidth',wline);
H=plot(t1,3*ones(size(t1)));
set(H,'Color',COL1(4,:),'LineWidth',wline);
set(gca,'Color','none','XColor','none','YColor','none','Xlim',t1([1 end]),'Ylim',[-.1 3.75])

axes(h11(2));
h11(2).clo;
H=plot(t1,zeros(size(t1)));
set(H,'Color',COL1(1,:),'LineWidth',wline);
hold on;
H=Vline(0,4,ctimes');
set(H,'Color',COL1(1,:),'LineWidth',2*wline,'LineStyle',':');
H=plot(t1,-stim+1);
set(H,'Color',COL1(2,:),'LineWidth',wline);
H=plot(t1,2*ones(size(t1)));
set(H,'Color',COL1(3,:),'LineWidth',wline);
H=plot(t1,3*ones(size(t1)));
set(H,'Color',COL1(4,:),'LineWidth',wline);
set(gca,'Color','none','XColor','none','YColor','none','Xlim',t1([1 end]),'Ylim',[-.1 3.75])

%probing
[h12,pos12] = tight_subplot(1,3,[0 dsmall],[0 0],[0 0],[],[],h1(3));
axes(h12(1));
h12(1).clo;
H=plot(t1,zeros(size(t1)));
set(H,'Color',COL1(1,:),'LineWidth',wline);
hold on;
H=Vline(0,4,ctimes');
set(H,'Color',COL1(1,:),'LineWidth',2*wline,'LineStyle',':');
H=plot(t1,zeros(size(stim))+1);
set(H,'Color',COL1(2,:),'LineWidth',wline);
H=plot(t1,remote+2);
set(H,'Color',COL1(3,:),'LineWidth',wline);
H=plot(t1,3*ones(size(t1)));
set(H,'Color',COL1(4,:),'LineWidth',wline);
set(gca,'Color','none','XColor','none','YColor','none','Xlim',t1([1 end]),'Ylim',[-.1 3.75])

axes(h12(2));
h12(2).clo;
H=plot(t1,zeros(size(t1)));
set(H,'Color',COL1(1,:),'LineWidth',wline);
hold on;
H=Vline(0,4,ctimes');
set(H,'Color',COL1(1,:),'LineWidth',2*wline,'LineStyle',':');
H=plot(t1,zeros(size(stim))+1);
set(H,'Color',COL1(2,:),'LineWidth',wline);
H=plot(t1,2*ones(size(t1)));
set(H,'Color',COL1(3,:),'LineWidth',wline);
H=plot(t1,3*ones(size(t1)));
set(H,'Color',COL1(4,:),'LineWidth',wline);
set(gca,'Color','none','XColor','none','YColor','none','Xlim',t1([1 end]),'Ylim',[-.1 3.75])

axes(h12(3));
h12(3).clo;
H=plot(t1,zeros(size(t1)));
set(H,'Color',COL1(1,:),'LineWidth',wline);
hold on;
H=Vline(0,4,ctimes');
set(H,'Color',COL1(1,:),'LineWidth',2*wline,'LineStyle',':');
H=plot(t1,zeros(size(stim))+1);
set(H,'Color',COL1(2,:),'LineWidth',wline);
H=plot(t1,2*ones(size(t1)));
set(H,'Color',COL1(3,:),'LineWidth',wline);
H=plot(t1,remote+3);
set(H,'Color',COL1(4,:),'LineWidth',wline);
set(gca,'Color','none','XColor','none','YColor','none','Xlim',t1([1 end]),'Ylim',[-.1 3.75])


%% second row - example rasters/PSTH
t=database.chin.t*1e3;
removebl=1;
if(removebl)
    ylim1=[-25 27];
    ylim2=[-8 8];
else
    ylim1=[0 45];
    ylim2=[0 25];
end
[hb,posb] = tight_subplot(2,2,[dbig dbig*1.5],[0 0],[0 0],[.4 .6],[],ha(2));

%% E cells
% [h20,pos20] = tight_subplot(1,2,[0 dbig],[0 0],[0 0],[],[],hb(1));

h21=plot_rast_psth({database.EcontRasters.pairA_rast,database.EcontRasters.probeA_rast,...
    database.EcontRasters.probeB_rast,database.EcontRasters.pairB_rast},hb(1),COL(3,:));
axes(h21(1))
plot_stims(1,1,0,COL1);
axes(h21(3))
plot_stims(0,1,0,COL1);
axes(h21(5))
plot_stims(0,0,0,COL1);
axes(h21(7));
plot_stims(1,0,0,COL1);
set(h21(7:8),'XColor',1-bgc);
set(h21([4 6]),'Ylim',[0 30],'YTick',[10 20]);
% % h21=plot_rast_psth({database.EcontRasters.pairA_rast,database.EcontRasters.pairB_rast,...
% %     database.EcontRasters.probeB_rast,...
% %     database.EcontRasters.probeA_rast},hb(1),COL(3,:));
% % axes(h21(1))
% % plot_stims(1,1,0,COL1);
% % axes(h21(3))
% % plot_stims(1,0,0,COL1);
% % axes(h21(5))
% % plot_stims(0,0,0,COL1);
% % axes(h21(7));
% % plot_stims(0,1,0,COL1);
% % set(h21(7:8),'XColor',1-bgc);
% % set(h21([6 8]),'Ylim',[0 30],'YTick',[10 20]);
% %rasters
% [h21,pos21] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h2(1));
% plot_rast_psth({database.EcontRasters.pairA_rast},h21(1),COL(3,:));
% plot_stims(1,1,0,COL1);
% % plot_rast_psth({database.Erasters.pairB_rast},h21(2),COL(3,:));
% plot_stims(1,0,0,COL1);
% set(h21(:),'XTick',[]);
% 
% [h23,pos23] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h2(3));
% plot_rast_psth({database.Erasters.probeA_rast},h23(1),COL(3,:));
% plot_stims(0,1,0,COL1);
% plot_rast_psth({database.Erasters.probeB_rast},h23(2),COL(3,:));
% set(h23(:),'XTick',[]);
% 
% [h25,pos25] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h2(5));
% plot_rast_psth({database.ECrasters.probeA_rast},h25(1),COL(3,:));
% plot_stims(0,0,1,COL1);
% plot_rast_psth({database.ECrasters.probeB_rast},h25(2),COL(3,:));
% set(h25(:),'XColor',[0 0 0],'TickLength',[.025 .025],'TickDir','both');
% 
[h2,pos2] = tight_subplot(3,1,[dsmall dbig],[0 0],[0 0],[],[],hb(3));

%pPSTH
[h22,pos22] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h2(1));
axes(h22(1));
plot_psth(t,[database.chin.pairA_E database.trunk.pairA_E],removebl,COL(3,:),ylim1);
axes(h22(2));
plot_psth(t,[database.chin.pairB_E database.trunk.pairB_E],removebl,COL(3,:),ylim1);
set(gca,'YColor','none');
set(h22(:),'XTick',[]);

[h24,pos24] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h2(2));
axes(h24(1));
plot_psth(t,[database.chin.probeA_E database.trunk.probeA_E],removebl,COL(3,:),ylim2);
axes(h24(2));
plot_psth(t,[database.chin.probeB_E database.trunk.probeB_E],removebl,COL(3,:),ylim2);
set(gca,'YColor','none');
set(h24(:),'XTick',[]);

[h26,pos26] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h2(3));
axes(h26(1));
plot_psth(t,[database.chin.probeC_E database.trunk.probeC_E],removebl,COL(3,:),ylim2);
axes(h26(2));
plot_psth(t,[database.chin.probeB_E database.trunk.probeB_E],removebl,COL(3,:),ylim2);
set(gca,'YColor','none');

%% I cells
removebl=1;
if(removebl)
    ylim1=[-10 27];
    ylim2=[-7 12];
else
    ylim1=[0 45];
    ylim2=[0 25];
end
% [h30,pos30] = tight_subplot(1,2,[0 dbig],[0 0],[0 0],[],[],hb(2));

h31=plot_rast_psth({database.IcontRasters.pairA_rast,database.IcontRasters.probeA_rast,...
    database.IcontRasters.probeB_rast,database.IcontRasters.pairB_rast},hb(2),COL(4,:));
axes(h31(1));
plot_stims(1,1,0,COL1);
axes(h31(3));
plot_stims(0,1,0,COL1);
axes(h31(5));
plot_stims(0,0,0,COL1);
axes(h31(7));
plot_stims(1,0,0,COL1);
set(h31(7:8),'XColor',1-bgc);

set(h31([2 8]),'Ylim',[0 100],'YTick',[40 80]);
set(h31([4 6]),'Ylim',[0 60],'YTick',[20 40]);
% [h3,pos3] = tight_subplot(3,2,[dsmall dbig],[0 0],[0 0],[],[],hb(2));

%rasters
% [h31,pos31] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h3(1));
% plot_rast_psth({database.Irasters.pairA_rast},h31(1),COL(4,:));
% plot_stims(1,1,0,COL1);
% plot_rast_psth({database.Irasters.pairB_rast},h31(2),COL(4,:));
% set(h31(:),'XTick',[]);
% plot_stims(1,0,0,COL1);
% 
% [h33,pos33] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h3(3));
% plot_rast_psth({database.Irasters.probeA_rast},h33(1),COL(4,:));
% set(gca,'Ylim',[0 80])
% plot_stims(0,1,0,COL1);
% plot_rast_psth({database.Irasters.probeB_rast},h33(2),COL(4,:));
% set(h33(:),'XTick',[]);
% set(gca,'Ylim',[0 80])
% 
% [h35,pos35] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h3(5));
% plot_rast_psth({database.ICrasters.probeA_rast},h35(1),COL(4,:));
% plot_stims(0,0,1,COL1);
% plot_rast_psth({database.ICrasters.probeB_rast},h35(2),COL(4,:));
% set(h35(:),'XColor',[0 0 0],'TickLength',[.025 .025],'TickDir','both');
% 
[h3,pos3] = tight_subplot(3,1,[dsmall dbig],[0 0],[0 0],[],[],hb(4));
%pPSTH
[h32,pos32] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h3(1));
axes(h32(1));
plot_psth(t,[database.chin.pairA_I database.trunk.pairA_I],removebl,COL(4,:),ylim1);
axes(h32(2));
plot_psth(t,[database.chin.pairB_I database.trunk.pairB_I],removebl,COL(4,:),ylim1);
set(gca,'YColor','none');
set(h32(:),'XTick',[]);

[h34,pos34] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h3(2));
axes(h34(1));
plot_psth(t,[database.chin.probeA_I database.trunk.probeA_I],removebl,COL(4,:),ylim2);
axes(h34(2));
plot_psth(t,[database.chin.probeB_I database.trunk.probeB_I],removebl,COL(4,:),ylim2);
set(gca,'YColor','none');
set(h34(:),'XTick',[]);

[h36,pos36] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h3(3));
axes(h36(1));
plot_psth(t,[database.chin.probeC_I database.trunk.probeC_I],removebl,COL(4,:),ylim2);
axes(h36(2));
plot_psth(t,[database.chin.probeB_I database.trunk.probeB_I],removebl,COL(4,:),ylim2);
set(gca,'YColor','none');

%%
return;
%%
%pairing
[h21,pos21] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h2(1));
plot_rast_psth({database.Irasters.pairA_rast database.Irasters.pairA_rast},h21(1),COL([3 4],:));
plot_rast_psth({database.Irasters.pairB_rast database.Irasters.pairB_rast},h21(2),COL([3 4],:));

%probing
[h22,pos22] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h2(2));
plot_rast_psth({database.Irasters.probeA_rast database.Irasters.probeA_rast},h22(1),COL([3 4],:));
plot_rast_psth({database.Irasters.probeB_rast database.Irasters.probeB_rast},h22(2),COL([3 4],:));


%% third row - population PSTHs
t=database.chin.t*1e3;
removebl=1;
if(removebl)
    ylim1=[-25 27];
    ylim2=[-10 10];
else
    ylim1=[0 45];
    ylim2=[0 25];
end
[h3,pos3] = tight_subplot(1,2,[0 dbig],[0 0],[0 0],[],[],ha(3));

%pairing
[h31,pos31] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h3(1));
axes(h31(1));
h31(1).clo;
plot_psth(t,[database.chin.pairA_E database.trunk.pairA_E],removebl,COL(3,:),ylim1);
plot_psth(t,[database.chin.pairA_I database.trunk.pairA_I],removebl,COL(4,:),ylim1);

axes(h31(2));
h31(2).clo;
plot_psth(t,[database.chin.pairB_E database.trunk.pairB_E],removebl,COL(3,:),ylim1);
plot_psth(t,[database.chin.pairB_I database.trunk.pairB_I],removebl,COL(4,:),ylim1);
set(gca,'YColor','none');

%probing
[h32,pos32] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h3(2));

axes(h32(1));
h32(1).clo;
plot_psth(t,[database.chin.probeA_E database.trunk.probeA_E],removebl,COL(3,:),ylim2);
plot_psth(t,[database.chin.probeA_I database.trunk.probeA_I],removebl,COL(4,:),ylim2);

axes(h32(2));
h32(2).clo;
plot_psth(t,[database.chin.probeB_E database.trunk.probeB_E],removebl,COL(3,:),ylim2);
plot_psth(t,[database.chin.probeB_I database.trunk.probeB_I],removebl,COL(4,:),ylim2);
set(gca,'YColor','none');

%controls
% [h33,pos33] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h3(3));

% axes(h33(1));
% h33(1).clo;
% plot_psth(t,database.ctl.probeA_E,removebl,COL(3,:),ylim2);
% plot_psth(t,database.ctl.probeA_I,removebl,COL(2,:),ylim2);
% 
% axes(h33(2));
% h33(2).clo;
% plot_psth(t,database.ctl.probeB_E,removebl,COL(3,:),ylim2);
% plot_psth(t,database.ctl.probeB_I,removebl,COL(2,:),ylim2);
% set(gca,'YColor','none');

% axes(h33(1));
% h33(1).clo;
% plot_psth(t,[database.chin.probeC_E database.trunk.probeC_E],removebl,COL(3,:),ylim2);
% plot_psth(t,[database.chin.probeC_I database.trunk.probeC_I],removebl,COL(2,:),ylim2);

% axes(h3(5));
% h3(5).clo;
% plot_psth(t,[database.chin.probeC_E database.trunk.probeC_E],removebl,COL(3,:),ylim2);
% plot_psth(t,[database.chin.probeC_I database.trunk.probeC_I],removebl,COL(2,:),ylim2);
% 
% axes(h3(6));
% h3(6).clo;
% plot_psth(t,[database.chin.indC_E database.trunk.indC_E],removebl,COL(3,:),ylim2);
% plot_psth(t,[database.chin.indC_I database.trunk.indC_I],removebl,COL(2,:),ylim2);

function plot_psth(t,Y,removebl,colvec,ylim)
fontsize=7;
nline=0.25;
wline=.5;
msize=3;

indbl=find(t<=10);
% Y=[database.chin.pairA_E database.trunk.pairA_E];
% Y=[database.nodelay.pairA_E]% database.trunk.pairA_I];
if(removebl)
    BL=mean(Y(indbl,:),1); %baseline firing
    Y=Y-ones(size(Y,1),1)*BL;
end

M=mean(Y,2);
SEM=std(Y,[],2)/sqrt(size(Y,2));
my_plotWithConf(t,M,SEM,colvec);
set(gca,'FontSize',fontsize,'XLim',t([1 end]),'Ylim',ylim,...
    'XTick',[0 50 100 150],'XTickLabelMode','auto','YTickLabelMode','auto','Color','none','TickLength',[.025 .025],'TickDir','both');
end


function strct=get_rasters(fname,modeC)
    S=load(fname);
    flds=fieldnames(S);
    cmd_field=flds{find(cellfun(@(x) numel(strfind(x,'Ch60')),flds))};
    spk_field=flds{find(cellfun(@(x) numel(strfind(x,'Ch9')),flds))};
    stim_field=flds{find(cellfun(@(x) numel(strfind(x,'Ch10')),flds))};
    if(modeC)
        cntxt_field=flds{find(cellfun(@(x) numel(strfind(x,'Ch7')),flds))};    
    else
        cntxt_field=flds{find(cellfun(@(x) numel(strfind(x,'Ch6')),flds))};    
    end
    
    cmd=S.(cmd_field).times;
    spk=S.(spk_field).times;
    cntxt=S.(cntxt_field).times;
    stim=S.(stim_field).times;    
    
    %align raster
    wind=[-0.020 .150];    
    clear strct;
    strct.pairA_rast=[]; pairAk=1;
    strct.pairB_rast=[]; pairBk=1;
    strct.probeA_rast=[];probeAk=1;    
    strct.probeB_rast=[];probeBk=1;
    
    for i=1:numel(cmd)
        pair=(sum(inrange(stim-cmd(i),wind))>0);
        block=sum(inrange(cntxt-cmd(i),wind));        
        ind_spk=find(inrange(spk-cmd(i),wind));
        if(numel(ind_spk))
            sptimes=spk(ind_spk)-cmd(i);
        end
        switch pair*10+block;
            case 0 %B probe
                strct.probeB_rast=[strct.probeB_rast;sptimes probeBk*ones(size(sptimes))];
                probeBk=probeBk+1;
            case 1 %A probe
                strct.probeA_rast=[strct.probeA_rast;sptimes probeAk*ones(size(sptimes))];
                probeAk=probeAk+1;
            case 10 %B pair
                strct.pairB_rast=[strct.pairB_rast;sptimes pairBk*ones(size(sptimes))];
                pairBk=pairBk+1;
            case 11 %A pair
                strct.pairA_rast=[strct.pairA_rast;sptimes pairAk*ones(size(sptimes))];
                pairAk=pairAk+1;
        end
    end     
end

function h0=plot_rast_psth(rasters,Ax,cols)
    fontsize=7;
    nline=0.25;
    wline=.5;
    msize=4;
    alpha=0.1; %smoothing kernel
    t0=[-19:2:143];
    t=[-5:2:143];
    bsize=2e-3;
    indt=find(t0>=t(1));
    edges=bin2edge(t0);
    if(size(cols,1)==1)
        cols=ones(numel(rasters),1)*cols;
    end
    
    %generate PSTH
    for r=1:numel(rasters)
        raster=rasters{r};
        K=max(raster(:,2));    
        h=zeros(K,numel(t0));
        for k=1:K
            spt=raster(raster(:,2)==k,1)*1e3;
            h0=histcounts(spt,edges);
            h(k,:)=filter(alpha,[1 alpha-1],h0);
        end        
        psths{r}=h(:,indt)'/bsize;
    end    
    
    [h0,pos0] = tight_subplot(numel(rasters),2,[0 0.075],[0 0],[0 0],[],[],Ax);
    
%     [h01,pos01] = tight_subplot(numel(rasters),1,[0 0],[0 0],[0 0],[],[],h0(1));
%     axes(Ax);
    for r=1:numel(rasters)
        raster=rasters{r};
        axes(h0(2*r-1));
        scatter(raster(:,1)*1e3,raster(:,2),msize,cols(r,:),'filled');
        set(gca,'Color','none','Xlim',t([1 end]),'FontSize',fontsize,'YTick',[],'YColor','none','XColor','none','Ylim',[0 80],'YDir','reverse');%max(raster(:,2))]);

        axes(h0(2*r));
        plot_psth(t,psths{r},0,cols(r,:),[0 50]); 
        set(gca,'Color','none','Xlim',t([1 end]),'FontSize',fontsize,...
            'YTick',[20 40],'XTick',[0 50 100],'XTickLabel',[]);
%         hold on;
    end
    set(gca,'XTickLabelMode','auto');
end

function plot_stims(stim,cntxt,ctl,cols)
    ylim=get(gca,'Ylim');
    W=3;
    hold on;
    if(stim)
        A=area(15+[0 W],ylim(2)*[1 1],ylim(1));
        set(A,'FaceColor',cols(2,:),'LineStyle','none');
    end
    if(cntxt)
        A=area(4.5+[0 W],ylim(2)*[1 1],ylim(1));
        set(A,'FaceColor',cols(3,:),'LineStyle','none');
    end
    if(ctl)
        A=area(4.5+[0 W],ylim(2)*[1 1],ylim(1));
        set(A,'FaceColor',cols(4,:),'LineStyle','none');
    end
    
        
end