D=load('model_delayed_GC3000.mat');
ex_model=D.strct1.ex_model;
models=D.strct1.models;
nettype=D.strct1.nettype;
nhidden=D.strct1.nhidden;
grid_M=D.strct1.grid_M;        
R=D.strct1.R;
TH=D.strct1.TH;

%% plot 
% ex_model=ex_model([1 2 4 5]);%for presentation
% models=models(:,[1 2 4 5]);

nline=0.25;
wline=0.5;
bwidth=4;
msize=.5;
fontsize=7;
N=size(models,1);
K=size(ex_model,2);
bgcol=[1 1 1];

Fall=figure;
BCOL=brewermap(12,'Paired');
COL=[[226 139 138]/255;.1216 .4706 .7059;[151 93 166;216 84 39;127 63 152]/255];

set(Fall,'Units','centimeters','OuterPosition',[0 1 17.2 10]);
[h0,pos0]=tight_subplot(2,1,[0.05 0],[0 0],[0 0],[.5 .5],[]);%Nh, Nw, [gap_h gap_w], [lower upper], [left right]
%% first row
[h1,p1]=tight_subplot(1,4,[0 0.04],[0 0],[0 0],[],[.35 .275 .175 .2],h0(1));%Nh, Nw, [gap_h gap_w], [lower upper], [left right]
h1(1).Visible=0;

axes(h1(2));
h1(2).clo;
for k=1:K
    m=ex_model(k).out_Re;
    n=ex_model(k).N;
    
    m=m(n>0);
    n=n(n>0);
    n=n/max(n)*25;

    S=scatter(k+randn(size(m))*.05,m,n,m,'filled','MarkerEdgeColor',0*[1 1 1],'MarkerFaceAlpha',1,'LineWidth',.1);
    hold on;
end
set(gca,'Xlim',[.5 K+.5],'Ylim',[-.5 .5],'Clim',[-.5 .5],'XAxisLocation','origin','XTick',[]);
if(bgcol==[1 1 1])
    set(gca,'Colormap',flipud(brewermap(64,'RdBu')));
else
    set(gca,'Colormap',invert_map(flipud(brewermap(64,'RdBu'))));    
end
set(gca,'XColor',1-bgcol,'YColor',1-bgcol,'Color','none','FontSize',fontsize);        
%% stats
for i=1:N
    re_var_flat(i,:)=cell2mat({models(i,:).var_re_flat});
    re_flat(i,:)=cell2mat({models(i,:).total_re_flat});
end
re_flat=re_flat./(mean(re_flat(:,1)));
X=reshape(re_flat,[],1);
G=reshape(ones(N,1)*[1:K],[],1);

axes(h1(3));
h1(3).clo;
boxplot(X,G,'PlotStyle','compact','Colors',COL,'Symbol','.w');
set(gca,'YLim',[0 1.35],'YTickMode','auto','YTickLabelMode','auto','XTickLabel',[],'Xlim',[0 K]+.5);
set(gca,'TickLength',[.025 .025],'TickDir','both');
set(gca,'Color','none','XColor',1-bgcol,'YColor',1-bgcol,'FontSize',fontsize,'box','off');
set(gca,'Position',p1{3});
bx=findobj('Tag','Box');
set(bx,'LineWidth',bwidth);
ws=findobj('Tag','Whisker');
set(ws,'LineWidth',wline);
cro=findobj('Tag','MedianInner');
set(cro,'Marker','none');
cri=findobj('Tag','MedianOuter');
set(cri,'MarkerFaceColor',[1 1 1],'MarkerSize',3);
%% NI
clear yy xx nn
xx=ex_model(1).out_Re(:);
% yy=reshape(cell2mat({ex_model([2 4 5]).ni_Re}),[],3);
yy(:,1)=ex_model(2).ni_Re(:);
yy(:,2)=ex_model(4).ni_Re(:);
yy(:,3)=ex_model(end).ni_Re(:);
nn=ex_model(1).N(:);
% xx=xx(nn>0); nn=nn(nn>0);
xx(nn==0)=nan;
nn(nn==0)=nan;
% yy(:,1)=yy(nn>0,1); yy(:,2)=yy(nn>0,2); yy(:,3)=yy(nn>0,3);
nn=nn/max(nn)*5;

axes(h1(4));
h1(4).clo;
I=[2 4 5];
for i=1:3
    S=scatter(xx,yy(:,i),nn,COL(I(i),:),'filled','MarkerEdgeColor',0*[1 1 1],'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.25,'LineWidth',.1);
    set(S,'MarkerFaceColor',COL(I(i),:));
    hold on;
end
plot([-.5 .5],[.5 -.5],'Color',.5*[1 1 1],'LineStyle','--');
set(gca,'Xlim',[-.5 .5],'Ylim',[-.5 .5],'Clim',[-.5 .5],'XAxisLocation','origin','XTick',[-.5 .5],'YTick',[-.5 .5],'YTickLabel',[-.5 .5],'YAxisLocation','origin');
set(gca,'XColor',1-bgcol,'YColor',1-bgcol,'Color','none');        
%% 2nd row: maps
[h2,p2]=tight_subplot(1,5,[0 0.01],[0 0],[0 0.05],[],[],h0(2));%Nh, Nw, [gap_h gap_w], [lower upper], [left right]
for i=1:5
    axes(h2(i));
    h2(i).clo;
    zz=ex_model(i).out_Re;
    nn=ex_model(i).N;
    Mzz=[zz zz(:,1)];
    AL=[nn nn(:,1)];
    
    S=polarplot3d(Mzz,AL,'AngularRange',[pi/2 5*pi/2],'plottype','surfa','RadialRange',[min(R(:)) max(R(:))],'AxisLocation','off','MeshScale',.5*[1 1]);       
    view(0,90);   
    hold on;  
    params.fish_length=19;
    params.fish_width=1.5;%cm
    params.tail_angle=0;
    params.tail_p=0.55; %non-bending portion of body     
    params.msize=18;
    params.mpos='';
    params.mneg='';
    params.bgcol=bgcol;
    
    [x]=get_skin(params);
    plot_skin(x,params);
    set(gca,'CLim',[-1 1],'Alim',[0 1.25],'XDir','reverse');
    set(gca,'Color','none','XColor','none','YColor','none','XGrid','off','YGrid','off');    
    axis('image');
end    

set(h2(:),'CLim',[-.5 .5]);
if(bgcol==[1 1 1])
    set(h2(:),'Colormap',flipud(brewermap(64,'RdBu')));
else
    set(h2(:),'Colormap',invert_map(flipud(brewermap(64,'RdBu'))));    
end
axes(h2(end));
C=colorbar;
set(C,'Location','east','AxisLocation','out','LineWidth',nline,'TickLength',0.025,'TickDirection','both','FontSize',fontsize,'FontWeight','normal');
p=C.Position;
C.Position(1)=p(1)+.025
%%
hall=findobj();
srf=findobj('Type','Surface');
ptch=findobj('Type','Patch');
set(srf,'Visible','off');
set(ptch,'Visible','off');
% export to metafile
print('-clipboard','-dmeta');

for i=2:numel(hall)
    if(~strcmp(hall(i).Type,'root') & ~strcmp(hall(i).Type,'figure'))
        hall(i).Visible='off';
    end
end
set(srf,'Visible','on');
set(ptch,'Visible','on');
print('-clipboard','-dbitmap','-r720');
