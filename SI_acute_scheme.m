% plot acute data figure

%% figure
% bgc=[0 0 0];
bgc=[1 1 1];
fontsize=7;
nline=0.25;
wline=.5;
msize=3;

Fall=figure;
set(Fall,'Units','centimeters');
set(Fall,'OuterPosition', [2  2  4  4],'Color',bgc);
COL=brewermap(10,'Set1');
COL(3,:)=[74 133 34]/255;
COL(4,:)=[112 130 56]/255;
COL1=[13 114 186;226 139 138;128 64 152;236 177 32]/255;
dbig=.075;
dsmall=.025;

% [ha, pos] = tight_subplot(2, 1, [.05 .0],[.05 .05],[.05 .0],[.2 .8],[]); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
%% first row - schemes
pwidth=.03;
T=1;
ctimes=[0.2 0.5 0.75];
stimes=ctimes+.01;
stimes_a=ctimes+.1;
stimes_b=ctimes+.01;
rtimes=ctimes+0.01;

t1=linspace(0,T,1e3);
stim=zeros(size(t1));
stim_a=zeros(size(t1));
stim_b=zeros(size(t1));
remote=stim;

for i=1:numel(ctimes)
    stim_a(inrange(t1-stimes_a(i),[0 pwidth]))=.5;
    stim_b(inrange(t1-stimes_b(i),[0 pwidth]))=.5;
    remote(inrange(t1-rtimes(i),[0 pwidth]))=.5;
end

[h1,pos1] = tight_subplot(1,2,[0 2*dsmall],[0 0],[0 0.0],[],[.4 .6]);
%pairing
[h11,pos11] = tight_subplot(1,2,[0 dsmall],[0 0],[0 0],[],[],h1(1));
axes(h11(1));
h11(1).clo;
H=plot(t1,zeros(size(t1)));
set(H,'Color',COL1(1,:),'LineWidth',wline);
hold on;
H=Vline(0,4,ctimes');
set(H,'Color',COL1(1,:),'LineWidth',2*wline,'LineStyle',':');
H=plot(t1,stim_a+1);
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
H=plot(t1,stim_b+1);
set(H,'Color',COL1(2,:),'LineWidth',wline);
H=plot(t1,2*ones(size(t1)));
set(H,'Color',COL1(3,:),'LineWidth',wline);
H=plot(t1,3*ones(size(t1)));
set(H,'Color',COL1(4,:),'LineWidth',wline);
set(gca,'Color','none','XColor','none','YColor','none','Xlim',t1([1 end]),'Ylim',[-.1 3.75])

%probing
[h12,pos12] = tight_subplot(1,3,[0 dsmall],[0 0],[0 0],[],[],h1(2));
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
