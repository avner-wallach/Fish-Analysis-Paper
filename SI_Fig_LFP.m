% function SI_Fig_LFP()
path(path,'Z:\GitHub\Fish-Model');

adir1='20200123'; seg1=11; chnum1=22; phi1=0.6*pi;
adir2='20190624'; seg2=1; chnum2=18; phi2=0.25*pi;
adir3='20190131'; seg3=2; chnum3=18; phi3=0.25*pi;
% script for visualizing freely moving recordings
apath='Z:\analysis_data';
if(~exist('eod1'))
    load([apath,filesep,adir1,filesep,'output'],'ops','eod','file');
    eod1=eod;
    ops1=ops;
    file1=file;
end
if(~exist('eod2'))
    load([apath,filesep,adir2,filesep,'output'],'ops','eod','file');
    eod2=eod;
    ops2=ops;
    file2=file;
end
if(~exist('eod3'))
    load([apath,filesep,adir3,filesep,'output'],'ops','eod','file');
    eod3=eod;
    ops3=ops;
    file3=file;
end
setenv('SESSDATE','20200107');
setenv('DATAPATH','Z:\mormyrid_data');
setenv('FRAMESCALE',num2str(ops1.framescale));

load('Z:\mormyrid_data\fish_images\fish_silhouette.mat');
global str_image bgc;
str_image=IMGS.cdata;
setenv('PXLSIZE',num2str(ops1.pxlsize)); %pxl to mm conversion
%%
bgc=[1 1 1];
d=0.05;

Fall=figure;
set(Fall,'Units','centimeters');
set(Fall,'OuterPosition', [1  2  17.2  10],'Color',bgc);
[ha, pa] = tight_subplot(1 , 4, [d d],[0 d],[0 0],[],[.19 .27 .27 .27]); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in
%% first column
% [hb, pb] = tight_subplot(2 , 1, [d 0],[0 0],[0 0],[],[],ha(1)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in
ha(1).Visible='off';
% params.plot_potential=1;
% params.plot_field=0;
% params.plot_lfield=1;
% params.reflection=0;

% axes(hb(1));
% hb(1).clo;
% plot_map('tail_angle',.225*pi,'tank_radius',23,'wall_dist',5,'wall_angle',.95*pi,...
%     'fish_length',15,'grid_M',50,'r_max',20,'tail_p',0.475,...
%     'plot_potential',0,'plot_field',0,'plot_lfield',0,'reflection',1,'mpos','.','mneg','.','axes',hb(1));
% axis('image');
% set(gca,'Xlim',[-20 8],'YLim',[-25 10]);

% axes(hb(1));
% hb(1).clo;
% if(~exist('frames'))
%     [frames,ax_grid] = get_fish_image(0,13259,'headfix','off','track','features','trackname','DLC_resnet50_fish_trackerSep11shuffle1_1030000');
% end
% Im=findobj('Type','Image');
% x0=621; y0=606; R=600;
% [X,Y]=meshgrid(1:size(Im.CData,2),1:size(Im.CData,1));
% blnkind=find(((X-x0).^2 + (Y-y0).^2)>R^2);
% CD=frames.cdata(:,:,1);
% CD(blnkind)=256;
% Im.CData=repmat(CD,1,1,3);
% set(gca,'Box','off','XColor','none','YColor','none');

% figure;
% A=axes;
% set(A,'Position',[0 0 1 1]);
% plot_map('tail_angle',.225*pi,'tank_radius',23,'wall_dist',5,'wall_angle',.95*pi,...
%     'fish_length',15,'grid_M',500,'grid_center',[18 -3],'r_max',30,'tail_p',0.475,...
%     'plot_potential',1,'plot_field',1,'plot_lfield',0,'reflection',0,'mpos','','mneg','','axes',A);
% colormap(A,lighter(brewermap(64,'BrBG'),1));
% axis('image');
% F=getframe(A);
% 
% axes(hb(3));
% hb(3).clo;
% imshow(F.cdata);

%% second row

% E=eod1(seg1);
% F=file1(seg1);
% Ch=chnum1;
% Phi=phi1;
% Ax=h(1);
[h1,p1]=data_vs_model(eod3(seg3),file3(seg3),15,0.7*pi,ha(2));
% h1(1).CLim=[-2 2.75]; h1(4).CLim=[-0.025 0.022];
[h2,p2]=data_vs_model(eod2(seg2),file2(seg2),chnum2,phi2,ha(3));
% h2(1).CLim=[-1.25 2];
h2(2).ALim=[1 2];
% h2(4).CLim=[-0.03 0.041];
[h3,p3]=data_vs_model(eod3(seg3),file3(seg3),16,1.25*pi,ha(4));
h3(2).ALim=[1 2];
% h3(1).CLim=[-1.75 2];
% h3(3).CLim=[-0.05 0.075];
% h3(4).CLim=[-0.032 0.03];
%% export
set(Fall,'Color','none');
hall=findobj();
for i=2:numel(hall)
     hall(i).Visible='on';
 end

srf=findobj('Type','Surface');
img=findobj('Type','Image');
ptch=findobj('Type','Patch');
set(srf,'Visible','off');
set(img,'Visible','off');
set(ptch,'Visible','off');
% export to metafile
print('-clipboard','-dmeta');

for i=2:numel(hall)
    if(~strcmp(hall(i).Type,'root') & ~strcmp(hall(i).Type,'figure'))
        hall(i).Visible='off';
    end
end
set(srf,'Visible','on');
set(img,'Visible','on');
set(ptch,'Visible','on');
print('-clipboard','-dbitmap','-r720');


%%
function [ha,pa]=data_vs_model(E,F,Ch,Phi,Ax)
    global str_image bgc;    
    fontsize=7;  
    clim=[-1 1];
    [ha, pa] = tight_subplot(2, 2, [0.05 0],[0 0],[0 0],[],[.3 .7],Ax); %Nh, Nw, [gap_h gap_w], [lower upper], [left right],dist_h,dist_w,axes_in

%     axes(ha(2));
%     ha(2).clo;
%     plot_tuning_polar(E,F,Ch,'t_col',6,'mfunc','offset','clim',[-2 2],...
%         'image',str_image,'fontsize',fontsize,'bgcol',bgc,...
%         'upsamp',0.25,'r_nbins',15,'th_nbins',15);
%     colorbar('off');
%     set(ha(1),'Position',pa{1});
    axes(ha(2));
    ha(2).clo;
    [M1,N1,hh]=plot_tuning_polar(E,F,Ch,'t_col',6,'mfunc','slope','clim',clim,...
        'image',str_image,'fontsize',fontsize,'bgcol',bgc,...
        'upsamp',0.25,'r_nbins',15,'th_nbins',15);
    colorbar('off');
    set(ha(2),'Position',pa{2});
    
    axes(ha(1));
    ha(1).clo;
    m=M1(:);
    n=N1(:);
    n=n/max(n)*50+1;
    S=scatter(randn(size(m))*.1,m,n,m,'filled','MarkerEdgeColor',0.25*[1 1 1],'MarkerFaceAlpha',1,'LineWidth',.1);
    set(ha(1),'Xlim',[-1 1],'Ylim',clim,'YTick',[-1 -.5 0 .5 1],'Clim',clim,'XAxisLocation','origin','XTick',[],'FontSize',fontsize);
    set(ha(1),'Colormap',ha(2).Colormap,'Color','none');

    ftemp=figure; atmp=axes;
    [FF,Z,T,R,TH,N,P1,P2]=get_tail_wall([atmp ha(4)],Phi,1);
    close(ftemp);
    ha(4).CLim=clim;
%     colorbar('Location','east','AxisLocation','out','LineWidth',nline,'TickLength',0.025,'TickDirection','both');
    axes(ha(3));
    ha(3).clo;
    m=P1(:);
    n=N(:);
    n=n/max(n)*10+1;
    S=scatter(randn(size(m))*.1,m,n,m,'filled','MarkerEdgeColor',0.25*[1 1 1],'MarkerFaceAlpha',1,'LineWidth',.1);
    set(ha(3),'Xlim',[-1 1],'Ylim',clim,'Clim',clim,'XAxisLocation','origin','XTick',[],'YTick',[-1 -.5 0 .5 1],'FontSize',fontsize);%,'Position',p3{4});
    set(ha(3),'Colormap',ha(4).Colormap,'Color','none');

%     ha(3).clo;  ha(4).clo;
%     get_tail_wall(ha([3 4]),Phi);
%     set(ha(3:4),'Xlim',[-28 28],'Ylim',[-23 23]);%,'Clim',[-0.05 0.03]);
%     axis(ha(3:4),'image');
%     set(ha(3),'XColor','none','YColor','none','Color','none','XGrid','off','YGrid','off');
end

