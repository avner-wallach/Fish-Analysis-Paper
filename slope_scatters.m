
F=figure;
set(gcf,'Color',[0 0 0]);
A=axes;
[Mout,N0,al]=plot_tuning_polar(eod2(2),file2(2),20,'t_col',6,'clim',[-1 1],...
'image',str_image.cdata,'fontsize',fontsize,'mfunc','slope','ops',ops2,'bgcol',bgc,...
'upsamp',0.25,'r_nbins',20,'th_nbins',20,'poly',p,'rastind',2);

m=Mout(:);
n=N0(:);
n=n/max(n)*25+1;
F=figure;
set(gcf,'Color',[0 0 0]);
A=axes;
S=scatter(randn(size(m))*.1,m,n,m,'filled','MarkerEdgeColor',0.25*[1 1 1],'MarkerFaceAlpha',1,'LineWidth',.1);
set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'Clim',[-1 1],'XAxisLocation','origin','XTick',[],'Position',[.1728 .11 .35 .815]);
colormap(gca,COL);
set(gca,'XColor',1-bgcol,'YColor',1-bgcol,'Color',bgcol);
set(gcf,'Position',[564 749 80 164]);

[F,Z,T,R,TH,N,P1,P2]=get_tail_wall([],.75*pi,1);
m=P1(:);
n=N(:);
n=n/max(n)*5+1;
F=figure;
set(gcf,'Color',[0 0 0]);
A=axes;
S=scatter(randn(size(m))*.1,m,n,m,'filled','MarkerEdgeColor',0.25*[1 1 1],'MarkerFaceAlpha',1,'LineWidth',.1);
set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'Clim',[-1 1],'XAxisLocation','origin','XTick',[],'Position',[.1728 .11 .35 .815]);
colormap(gca,COL);
set(gca,'XColor',1-bgcol,'YColor',1-bgcol,'Color',bgcol);
set(gcf,'Position',[564 749 80 164]);

row{1}=[2;2;2;3;3;3;4]; col{1}=[9;10;11;9;10;11;11];
row{2}=[16;16;16;17;17;17;18;18;18]; col{2}=[14;15;16;14;15;16;14;15;16];
row{3}=[10;10;10;11;11;11;12;12;12]; col{3}=[9;10;11;9;10;11;9;10;11];
for i=1:3
    z{i}=[];
    for j=1:numel(row{i})
        z{i}(j,:)=permute(Z(col{i}(j),row{i}(j),:),[3 1 2]);
    end
end

F=figure;
clim=[-.5 .5];
CC=brewermap(64,'Greys');    CCN=size(CC,1);
if(bgcol==[0 0 0])
    CC=flipud(CC);
end
d(1)=subplot(1,2,1);
d(2)=subplot(1,2,2);
[h1,p1]=tight_subplot(3, 1, [.05 .0],[.0 .00],[.00 .0],[],[],d(1)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
set(gcf,'Color',bgcol);
for i=1:3
    axes(h1(i));
    TT=repmat(T',size(z{i},1),1);
    S=scatter(TT(:),z{i}(:),10,TT(:),'filled');
    set(S,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',1);%,'MarkerFaceColor',COL(1,:));    
    hold on;    
    f=fit(TT(~isnan(z{i})),z{i}(~isnan(z{i})),'poly1');
    sl1(i)=f.p1;
    cval=min(max((f.p1-clim(1))/diff(clim),0),1);
    MCOL(i,:)=interp1(linspace(0,1,CN),COL,cval);
    XX=linspace(min(T),max(T),1e2);
    hold on;
    H=plot(XX,f(XX)); set(H,'Color',MCOL(i,:),'LineWidth',3); legend('off');
    
    legend('off');    
    xlabel('');
    ylabel('');
    set(gca,'XAxisLocation','origin','YAxisLocation','origin',...
        'Ylim',[-.5 .5],'Xlim',[-1.5 1.5],'YTick',[],'XTick',[]);    
    set(gca,'Color',bgcol,'YColor',1-bgcol,'XColor',1-bgcol,'Colormap',CC);    
end

axes(d(2));
for i=1:3
    H=plot(i,sl1(i),'o');
    set(H,'MarkerSize',10,'MarkerFaceColor',MCOL(i,:),'MarkerEdgeColor',[1 1 1],'LineWidth',1);
    hold on;
end
  set(gca,'Color',bgcol,'YColor',1-bgcol,'XColor',1-bgcol,'Xlim',[.5 3.5],'Ylim',[-.4 .5],'XTick',[],'YTick',[],'XAxisLocation','origin','Box','off');    
