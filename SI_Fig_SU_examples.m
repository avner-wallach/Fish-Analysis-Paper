adir1='20200123'; seg1=4;
adir2='20190624'; seg2=6;
adir3='20190131'; seg3=4;
adir4='20200709'; seg4=7;
adir5='20200113'; seg5=4; seg6=5;
fontsize=7;       

% adir6='20210817'; seg6=5;
apath='Z:\analysis_data';
if(~exist('eod1'))
    load([apath,filesep,adir1,filesep,'output'],'ops','eod','file');
    eod1=eod(seg1);
    ops1=ops;
    file1=file(seg1);
end
if(~exist('eod2'))
    load([apath,filesep,adir2,filesep,'output'],'ops','eod','file');
    eod2=eod(seg2);
    ops2=ops;
    file2=file(seg2);
end
if(~exist('eod3'))
    load([apath,filesep,adir3,filesep,'output'],'ops','eod','file');
    eod3=eod(seg3);
    ops3=ops;
    file3=file(seg3);
end
if(~exist('eod4'))
    load([apath,filesep,adir4,filesep,'output'],'ops','eod','file');
    eod4=eod(seg4);
    ops4=ops;
    file4=file(seg4);
end
if(~exist('eod5'))
    load([apath,filesep,adir5,filesep,'output'],'ops','eod','file');
    eod5=eod(seg5);
    ops5=ops;
    file5=file(seg5);
    eod6=eod(seg6);
    ops6=ops;
    file6=file(seg6);    
end
% if(~exist('eod6'))
%     load([apath,filesep,adir5,filesep,'output'],'ops','eod','file');
%     eod6=eod(seg6);
%     ops6=ops;
%     file6=file(seg6);
% end

load('Z:\mormyrid_data\fish_images\fish_silhouette.mat');
global str_image
str_image=IMGS.cdata;
setenv('PXLSIZE',num2str(ops1.pxlsize)); %pxl to mm conversion

% bgc=[0 0 0];
d=0.025;
COL=[236 157 118;210 149 0;74 133 34]/255;

Fall=figure;
set(Fall,'Units','centimeters');
set(Fall,'OuterPosition', [1  2  16.5  27],'Color',[1 1 1]);
[h, pos] = tight_subplot(7 , 1, [d d],[d 2*d],[2*d 2*d]); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in
%type 1
[h1,p1]=plot_unit(eod1,file1,23,ops1.seg(seg1).ind_lim(1,:),19,eod1.raster{1},COL(2,:),h(1),25,25,[10 70],'2d');
for i=[4 6]
    axes(h1(i));
    cb(i)=colorbar('northoutside','FontSize',fontsize);
    h1(i).Position=p1{i};
    cbp{i}=cb(i).Position;
    cb(i).Position=cbp{i}+[d 0 -2*d 0];    
    cb(i).TickDirection="both";
    cb(i).Ruler.TickLabelGapOffset = -1
end
[h2,p2]=plot_unit(eod5,file5,21,ops5.seg(seg5).ind_lim(3,:),17,eod5.raster{3},COL(2,:),h(2),25,25,[20 110],'2d');
% [h2,p2]=plot_unit(eod6,file6,24,ops6.seg(seg6).ind_lim(2,:),18,h(2));
% [h2,p2]=plot_unit(eod3,file3,19,ops3.seg(seg3).ind_lim(1,:),17,h(2));
% h2(1).CLim=[0 .65];
[h3,p3]=plot_unit(eod4,file4,23,ops4.seg(seg4).ind_lim(1,:),17,eod4.raster{1},COL(2,:),h(3),25,25,[0 40],'2d');
% h3(1).CLim=[0 1.1];
%type 2
[h4,p4]=plot_unit(eod4,file4,24,ops4.seg(seg4).ind_lim(2,:),17,eod4.raster{2},COL(3,:),h(4),25,25,[5 20],'2d');
% h5(1).CLim=[0 .75];
% [h6,p6]=plot_unit(eod5,file5,20,ops5.seg(seg5).ind_lim(2,:),17,h(6));
[h5,p5]=plot_unit(eod6,file6,22,ops6.seg(seg6).ind_lim(4,:),18,eod6.raster{4},COL(3,:),h(5),30,35,[2 11],'2d');
[h6,p6]=plot_unit(eod2,file2,23,ops2.seg(seg2).ind_lim(1,:),21,eod2.raster{1},COL(3,:),h(6),25,25,[2 11],'2d');
[h7,p7]=plot_unit(eod2,file2,23,ops2.seg(seg2).ind_lim(1,:),21,eod2.raster{1},COL(3,:),h(7),10,25,[2 12],'pl');
% [h7,p7]=plot_unit_polar(eod2,file2,23,ops2.seg(seg2).ind_lim(1,:),21,eod2.raster{1},COL(3,:),h(7),10,25);
set(gca,'Xtick',[-10 0 10 20],'XTickLabel',[-10 0 10 20]);
%% export
set(Fall,'Color','none');
hall=findobj();
% for i=2:numel(hall)
%     hall(i).Visible='on';
% end

srf=findobj('Type','Surface');
img=findobj('Type','Image');
% rct=findobj('Type','Rectangle')
set(srf,'Visible','off');
set(img,'Visible','off');
% set(rct,'Visible','off');
% export to metafile
print('-clipboard','-dmeta');

for i=2:numel(hall)
    if(~strcmp(hall(i).Type,'root') & ~strcmp(hall(i).Type,'figure'))
        hall(i).Visible='off';
    end
end
set(srf,'Visible','on');
set(img,'Visible','on');
% set(rct,'Visible','on');
print('-clipboard','-dbitmap','-r720');


%%
function [hout,pout]=plot_unit(E,F,Ch,Ind,Lfp,Rast,C,Ax,Nx,Ny,clim,mode)
    global str_image;
    bgcol=[1 1 1];    
    fontsize=7;       
    msize=1;
    wline=1;
    Cslope=[-.5 .5];    
    xlim=[-550 550];
    ylim=[-550 450];
    r_nbins=10;
    th_nbins=25;
    xnbins=Nx;
    ynbins=Ny;
    rmax=500;
    upsamp=.25;
    minsamp=50;
   
    [ha, pa] = tight_subplot(1, 4, [0 0.045],[0 0],[0 0],1,[.13 .29 .29 .29],Ax); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in

 
    [h1, p1] = tight_subplot(1, 2, [0 0],[0 0],[0 0],1,[.2 .8],ha(2)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in
    axes(h1(2));
    [M1,N1,hh]=plot_tuning_polar(E,F,Ch,'mfunc','slope','t_col',6,'image',str_image,...
        'clim',Cslope/2,'bgcol',bgcol,'r_nbins',r_nbins,'th_nbins',th_nbins,...
        'minsamp',minsamp,'upsamp',upsamp,'ind_lim',Ind);
    colorbar('off');
    set(h1(2),'Alim',[1 2],'Color','none');

    axes(h1(1));
    h1(1).clo;
    m=M1(:);
    n=N1(:);
    n(n<1000)=nan;
    n=n/max(n)*50+1;
    S=scatter(randn(size(m))*.1,m,n,m,'filled','MarkerEdgeColor',0.25*[1 1 1],'MarkerFaceAlpha',1,'LineWidth',.1);
    hold on;
    mas=nansum(abs(M1(:)).*N1(:))/nansum(N1(:));
    T=text(0,.2,[num2str(mas,'%0.2f')],'FontSize',fontsize,'Color',C);
    set(h1(1),'Xlim',[-1 1],'Ylim',Cslope/2,'YTick',[-.25 0 .25],'Clim',Cslope/2,'XAxisLocation','origin','XTick',[],'FontSize',fontsize);
    set(h1(1),'Colormap',h1(2).Colormap,'Color','none');
    
    [h2, p2] = tight_subplot(1, 2, [0 0],[0 0],[0 0],1,[.2 .8],ha(3)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in
    axes(h2(2));
    [M2,N2,hh]=plot_tuning_polar(E,F,Ch,'mfunc','slope','t_col',Lfp,'image',str_image,...
        'clim',Cslope,'bgcol',bgcol,'r_nbins',r_nbins,'th_nbins',th_nbins,...
        'minsamp',minsamp,'upsamp',upsamp,'ind_lim',Ind,'r_max',rmax);
    colorbar('off');
    set(h2(2),'Alim',[1 2],'Color','none');
    
    axes(h2(1));
    h2(1).clo;
    m=M2(:);
    n=N2(:); 
    n(n<1000)=nan;
    n=n/max(n)*50+1;
    S=scatter(randn(size(m))*.1,m,n,m,'filled','MarkerEdgeColor',0.25*[1 1 1],'MarkerFaceAlpha',1,'LineWidth',.1);
    hold on;
    mas=nansum(abs(M2(:)).*N2(:))/nansum(N2(:));
    T=text(0,.45,[num2str(mas,'%0.2f')],'FontSize',fontsize,'Color',C);

    set(h2(1),'Xlim',[-1 1],'Ylim',Cslope,'YTick',[-1 -.5 0 .5 1],'Clim',Cslope,'XAxisLocation','origin','XTick',[],'FontSize',fontsize);
    set(h2(1),'Colormap',h2(2).Colormap,'Color','none');
    
    [h3, p3] = tight_subplot(1, 2, [0 0],[0 0],[0 0],1,[.2 .8],ha(4)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in    
    axes(h3(2));    
    if(mode=='2d')
        [M3,N3,hh]=plot_tuning2d(E,F,Ch,'mfunc','mean','t_col',6,'image',str_image,...
            'clim',clim,'bgcol',bgcol,'upsamp',1.25,'ind_lim',Ind,'zscale',[0 40],...
            'x_lim',xlim,'y_lim',ylim,'minsamp',minsamp,'r_max',rmax,'x_nbins',xnbins,'y_nbins',ynbins);
    else
        [M3,N3,hh]=plot_tuning_polar(E,F,Ch,'mfunc','mean','t_col',6,'image',str_image,...
            'clim',clim,'bgcol',bgcol,'upsamp',upsamp,'ind_lim',Ind,'zscale',[0 40],...
            'x_lim',xlim,'y_lim',ylim,'minsamp',minsamp,'r_max',rmax,'r_nbins',r_nbins,'th_nbins',th_nbins);
    end        
    colorbar('off');
    [cb]=colorbar('east','FontSize',fontsize,'AxisLocation','out','TickLength',0.025,'TickDirection','both');
    pb=cb.Position;
    cb.Position(1)=pb(1)+0.025;
    
    BW=1;
    m3=discretize(M3(~isnan(M3)),50);
    f=ksdensity(m3,[1:50],'Bandwidth',BW);
    [h,s]=entropy(m3);
    x=linspace(0,1,50);
    axes(h3(1));
    h3(1).clo;
    fp=4;
    P=patch([-f*fp fliplr(f*fp)],[x fliplr(x)],C);
    text(-.5,1,['S=',num2str(s,'%0.2f')],'FontSize',fontsize,'Color',C);
    set(P,'FaceColor',C,'LineStyle','none');
    set(h3(1),'FontSize',fontsize,'Color','none','Box','off','Xlim',[-.5 .5],'XTick',[],'YTickLabelMode','auto');

    bpre=-.5;
    bpost=1;
    rpre=10;
    rpost=25;
    binsize=1;
    ksmooth=3;

   %ha(4).Position = pa{4} + [0.02 0 -0.02 0]; %shift left
    [h0,p0] = tight_subplot(2, 1, [.00 0],[.0 .0],[.0 .025],[.35 .65],[1],ha(1)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in
    N=1e3;
    N0=size(E.data,1);
    axes(h0(1));
    h0(1).clo;    
    Rast(:,1)=Rast(:,1)*1e3;
%     r=randi(N0,N,1);
    r=find(E.data(:,1)>0.05);
    r=r(randperm(numel(r),N));    
    resp=nanmean(E.data(r,Ch))*40;    
    ind=find(ismember(Rast(:,2),r));        
    q=zeros(N0,1);
    q(r)=[1:N];
    drast=Rast(ind,:); 
    drast(:,2)=q(drast(:,2));
    rastplot(drast,'|',C,msize);  
    hold on;
    A=area([-bpre bpost],[N N],'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',0.2,'ShowBaseLine','off');
    set(h0(1),'Color',bgcol,'XColor','none','YColor',1-bgcol,'XTick',[],'YTick',[],'box','off');
    set(h0(1),'Xlim',[-10 25],'Ylim',[0 N]);    
    
    edges=[-rpre:binsize:(-bpre)  bpost:binsize:rpost];
    bins=edge2bin(edges);        
    indblank=find(inrange(bins,[-bpre bpost]));
    axes(h0(2));
    h0(2).clo;        
    h=histcounts(drast(:,1),edges)/N/binsize*1e3;

    hh=smooth(h,ksmooth);
    hh(indblank)=nan;
    H=plot(bins,hh);
    set(H,'Color',C,'LineWidth',wline);
    hold on;
    A=area([-bpre bpost],[1e3 1e3],'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',0.2);
    
    T=text(10,30,['R=',num2str(round(resp)),'Hz'],'FontSize',fontsize,'FontName','Arial','Color',C);    
    
    set(h0(2),'Color',bgcol,'XColor',1-bgcol,'YColor',1-bgcol,'FontSize',fontsize,'box','off','YAxisLocation','left');
    if(nanmax(hh)>60)
        set(gca,'Xlim',[-10 25],'Ylim',[0 250],'XTick',[-10 0 10 20],'XTickLabel',[],'YTick',[100 200]);        
    else
        set(gca,'Xlim',[-10 25],'Ylim',[0 60],'XTick',[-10 0 10 20],'XTickLabel',[],'YTick',[20 40]);                
    end
    
%     set(ha(:),'YLim',[-210 125]); 
    hout=[h0;h1;h2;h3];
    pout={p0{:} p1{:} p2{:} p3{:}}';    
    
end

function [ha,pa]=plot_unit_polar(E,F,Ch,Ind,Lfp,Rast,C,Ax,Nr,Nth)
    global str_image;
    bgcol=[1 1 1];    
    fontsize=7;       
    msize=1;
    wline=1;
    Cslope=[-.5 .5];    
    xlim=[-550 550];
    ylim=[-550 450];
    rnbins=Nr;
    thnbins=Nth;
    rmax=500;
    upsamp=.75;
    minsamp=20;
    tlength=0.025;

    K=4;
    [ha, pa] = tight_subplot(1, K, [0 0],[0 0],[0 0],1,ones(1,K)/K,Ax); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in

    axes(ha(1));
    plot_tuning_polar(E,F,Ch,'mfunc','offset','t_col',6,'image',str_image,...
        'clim',[-1 1],'bgcol',bgcol,'upsamp',upsamp,'ind_lim',Ind,...
        'minsamp',minsamp,'r_max',rmax,'r_nbins',rnbins,'th_nbins',thnbins,'objidx',0);
    colorbar('off');

    axes(ha(2));
    plot_tuning_polar(E,F,Ch,'mfunc','slope','t_col',6,'image',str_image,...
        'clim',Cslope/2,'bgcol',bgcol,'x_lim',xlim,'y_lim',ylim,...
        'minsamp',minsamp,'upsamp',upsamp,'ind_lim',Ind,'r_max',rmax,'objidx',0);
    colorbar('off');
%     axes(ha(3));
%     plot_tuning_polar(E,F,Ch,'mfunc','slope','t_col',1,'image',str_image,...
%         'clim',Cslope,'bgcol',bgcol,'x_lim',xlim,'y_lim',ylim,...
%         'minsamp',minsamp,'upsamp',upsamp,'ind_lim',Ind,'r_max',rmax,'objidx',0);
    colorbar('off');
    axes(ha(3));
    plot_tuning_polar(E,F,Ch,'mfunc','slope','t_col',Lfp,'image',str_image,...
        'clim',Cslope,'bgcol',bgcol,'x_lim',xlim,'y_lim',ylim,...
        'minsamp',minsamp,'upsamp',upsamp,'ind_lim',Ind,'r_max',rmax,'objidx',0);
    colorbar('off');

    bpre=-.5;
    bpost=1;
    rpre=10;
    rpost=25;
    binsize=1;
    ksmooth=3;

    ha(4).Position = pa{4} + [0.02 0 -0.02 0]; %shift left
    [hb1,pb1] = tight_subplot(2, 1, [.00 0],[.0 .0],[.0 .0],[.5 .5],[1],ha(4)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in
    N=1e3;
    N0=size(E.data,1);
    axes(hb1(1));
    hb1(1).clo;    
    Rast(:,1)=Rast(:,1)*1e3;
%     r=randi(N0,N,1);
    r=find(E.data(:,1)>0.05);
    r=r(randperm(numel(r),N));    
    resp=nanmean(E.data(r,Ch))*40;    

    ind=find(ismember(Rast(:,2),r));        
    
    q=zeros(N0,1);
    q(r)=[1:N];
    drast=Rast(ind,:); 
    drast(:,2)=q(drast(:,2));
    rastplot(drast,'|',C,msize);  
    hold on;
    A=area([-bpre bpost],[N N],'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',0.2);
    set(hb1(1),'Color',bgcol,'XColor',1-bgcol,'YColor',1-bgcol,'XTick',[],'YTick',[],'box','off');
    set(hb1(1),'Xlim',[-10 25],'Ylim',[0 N]);
    set(gca,'TickLength',[tlength tlength],'TickDir','both');

    
    edges=[-rpre:binsize:(-bpre)  bpost:binsize:rpost];
    bins=edge2bin(edges);        
    indblank=find(inrange(bins,[-bpre bpost]));
    axes(hb1(2));
    hb1(2).clo;        
    h=histcounts(drast(:,1),edges)/N/binsize*1e3;

    hh=smooth(h,ksmooth);
    hh(indblank)=nan;
    H=plot(bins,hh);
    set(H,'Color',C,'LineWidth',wline);
    hold on;
    A=area([-bpre bpost],[1e3 1e3],'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',0.2);
    
    T=text(0,10,[num2str(round(resp)),'Hz'],'FontSize',5,'FontName','Arial','Color',C);

    set(hb1(2),'Color','none','XColor',1-bgcol,'YColor',1-bgcol,'FontSize',fontsize,'box','off','YAxisLocation','left');
    if(nanmax(hh)>60)
        set(gca,'Xlim',[-10 25],'Ylim',[0 250],'XTick',[],'YTick',[100 200]);        
    else
        set(gca,'Xlim',[-10 25],'Ylim',[0 60],'XTick',[],'YTick',[20 40]);                
    end
    set(gca,'TickLength',[tlength tlength],'TickDir','both');

%     set(ha(:),'YLim',[-210 125]); 
end

