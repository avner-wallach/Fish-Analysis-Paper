% make panels for Figure 1: Afferent input in freely moving fish
%% load data
path(path,'Z:\GitHub\Fish-Model');
path(path,'Z:\GitHub\Fish-Pipeline')
if(~exist('eod'))
    load('Z:\analysis_data\20190131\output.mat')
end
load('Z:\mormyrid_data\fish_images\straight_image.mat');
str_image=IMGS;
setenv('PXLSIZE','0.423'); %pxl to mm conversion
pxl2mm=0.423;
setenv('FRAMESCALE',num2str(ops.framescale));
setenv('SESSDATE','20190131');
setenv('DATAPATH','Z:\mormyrid_data');
%% get snapshot images
getimages=0;
if(getimages)
    [r,th]=convert_exo_polar(frame(2).data(:,1),frame(2).data(:,2),...
        frame(2).data(:,3),file(2).circle,pxl2mm);
    % ind=find(inrange(r,[49 51]) & inrange(th,pi/2+[-pi/20 pi/20]));
    x=r.*sin(th); y=r.*cos(th);
    t=frame(2).data(:,9);
    for i=2:3
        figure(Fre);
        h=impoly();
        p=h.getPosition;
        ind = find(inpolygon(x,y,p(:,1),p(:,2)));
        %get indices for frames
        [fileidx{i},frameidx{i}]=get_frame_indices(ind,frame(2).ind0);
        T=mat2cell(t(ind),cellfun(@numel,frameidx{i}));
        [m,im]=max(cellfun(@numel,frameidx{i}));
        figure;
        [frames{i},ax_grid] = get_fish_image(ops.seg(2).files(fileidx{i}(im)),frameidx{i}{im},'headfix','off');
    end
    save('fish_frames','frames','fileidx','frameidx','-v7.3');
else
    load('fish_frames','frames','fileidx','frameidx','p');
end
%% start figure
Fall=figure;
bgc=[0 0 0];
% bgc=[1 1 1];
fontsize=16;
[ha, pos] = tight_subplot(3, 3, [.05 .05],[.05 .05],[.05 .05]);
% ha(1).Visible='off'; %place for moodel scheme

%% get joint figure of all fish
F(:,:,:,1)=frames{1}(10).cdata;
F(:,:,:,2)=frames{2}(10).cdata;
F(:,:,:,3)=frames{3}(14).cdata;
D=double(F);
th=45;
M=median(D,4);
x=1:size(M,2); y=1:size(M,1);
[X,Y]=meshgrid(x,y);
x0=170; y0=167; r=size(M,1)/2;
B0=(((X-x0).^2+(Y-y0).^2)<r.^2);
B1=(abs(D(:,:,:,1)-M)>th);
B2=(abs(D(:,:,:,2)-M)>th);
B3=(abs(D(:,:,:,3)-M)>th);
c1(:,:,1)=57/255*ones(size(B1(:,:,1))); c1(:,:,2)=116/255*ones(size(B1(:,:,1))); c1(:,:,3)=254/255*ones(size(B1(:,:,1)));
c2(:,:,1)=64/255*ones(size(B1(:,:,1))); c3(:,:,2)=202/255*ones(size(B1(:,:,1))); c3(:,:,3)=141/255*ones(size(B1(:,:,1)));
c3(:,:,1)=203/255*ones(size(B1(:,:,1))); c2(:,:,2)=193/255*ones(size(B1(:,:,1))); c2(:,:,3)=41/255*ones(size(B1(:,:,1)));
FF1=(255-D(:,:,:,1)).*c1;
FF2=(255-D(:,:,:,2)).*c2;
FF3=(255-D(:,:,:,3)).*c3;
axes(ha(2));
I=image(uint8(M));
set(I,'AlphaData',B0);
hold on;
I=image(uint8(FF1));
set(I,'AlphaData',B1(:,:,1).*B0);
I=image(uint8(FF2));
set(I,'AlphaData',B2(:,:,1).*B0);
I=image(uint8(FF3));
set(I,'AlphaData',B3(:,:,1).*B0);
set(ha(2),'XColor','none','YColor','none');
axis('image');
%% model re- and ex- afference examples
ha(1).Position=pos{1}+[0 0 0 -0.05];
[hb,pb]=split_axes(ha(1),2,3,0.02,0.02);
plot_map('tail_angle',0,'tank_radius',23,'wall_dist',15,'wall_angle',1*pi,'mpos','','mneg','','axes',hb(1));
colormap(hb(1),lighter(brewermap(64,'BrBG'),0.4));
plot_map('tail_angle',0,'tank_radius',23,'wall_dist',15,'wall_angle',.5*pi,'mpos','','mneg','','axes',hb(2),'grid_center',[0 0]);
colormap(hb(2),lighter(brewermap(64,'BrBG'),0.4));
plot_map('tail_angle',0,'tank_radius',23,'wall_dist',5,'wall_angle',pi/4,'mpos','','mneg','','axes',hb(3));
colormap(hb(3),lighter(brewermap(64,'BrBG'),0.4));
plot_map('tail_angle',pi/4,'tank_radius',2e3,'wall_dist',200,'wall_angle',0,'mpos','','mneg','','axes',hb(4));
colormap(hb(4),lighter(brewermap(64,'BrBG'),0.4));
plot_map('tail_angle',0,'tank_radius',2e3,'wall_dist',200,'wall_angle',0,'mpos','','mneg','','axes',hb(5));
colormap(hb(5),lighter(brewermap(64,'BrBG'),0.4));
plot_map('tail_angle',-pi/4,'tank_radius',2e3,'wall_dist',200,'wall_angle',0,'mpos','','mneg','','axes',hb(6));
colormap(hb(6),lighter(brewermap(64,'BrBG'),0.4));
%% rois
% pl{1}=[ -46.7480  -18.6736 ; -47.3649  -13.9881 ; -49.1734   -9.6219 ;  -52.0504   -5.8726 ;  -55.7997   -2.9956 ;...
%         -60.1659   -1.1870 ; -64.8514   -0.5702 ;  -69.5369   -1.1870 ;  -73.9031   -2.9956 ;  -77.6524   -5.8726 ;...
%         -80.5294   -9.6219 ; -82.3379  -13.9881 ;  -82.9548  -18.6736 ;  -82.3379  -23.3591 ;  -80.5294  -27.7253 ;...
%         -77.6524  -31.4746 ; -73.9031  -34.3516 ;  -69.5369  -36.1601 ;  -64.8514  -36.7770 ;  -60.1659  -36.1601 ;...
%         -55.7997  -34.3516 ; -52.0504  -31.4746 ; -49.1734  -27.7253 ;  -47.3649  -23.3591 ;  -46.7480  -18.6736];
% pl{2}=[ 31.1653  122.0287 ;   30.5669  124.7732 ;   28.7131  127.4077 ;   25.5319  129.7008 ;   21.1521  131.3027 ;...    
%    16.0410  131.8825 ;   10.9298  131.3027 ;    6.5501  129.7008 ;    3.3689  127.4077 ;    1.5151  124.7732 ;...
%     0.9166  122.0287 ;    1.5151  119.2843 ;    3.3689  116.6497 ;    6.5501  114.3566 ;   10.9298  112.7547 ;...
%    16.0410  112.1750 ;   21.1521  112.7547 ;   25.5319  114.3566 ;   28.7131  116.6497 ;   30.5669  119.2843 ;...
%    31.1653  122.0287];
% pl{3}=[  159.0349  -18.9027 ;  158.5759  -12.5436 ;  157.3031   -6.9655 ;  155.4767   -2.6844 ;  153.3930    0.1567 ;...
%   151.2689    1.7010 ;  149.1811    2.1797 ;  147.0933    1.7010 ;  144.9692    0.1567 ;  142.8855   -2.6844 ;...
%   141.0591   -6.9655 ;  139.7863  -12.5436 ;  139.3274  -18.9027 ;  139.7863  -25.2619 ;  141.0591  -30.8399 ;...
%   142.8855  -35.1211 ;  144.9692  -37.9622 ;  147.0933  -39.5065 ;  149.1811  -39.9852 ;  151.2689  -39.5065 ;...
%   153.3930  -37.9622 ;  155.4767  -35.1211 ;  157.3031  -30.8399 ;  158.5759  -25.2619 ;  159.0349  -18.9027];
%% plot ex afference, re afference maps
axes(ha(5));
ha(5).clo;
plot_tuning_polar(eod(2),file(2),18,'t_col',6,'clim',[-2 2],...
    'image',IMGS.cdata,'mfunc','mean','ops',ops,'bgcol',bgc,...
    'upsamp',0.25,'r_nbins',15,'th_nbins',15);
ha(5).Position=pos{5}-[0.055 0 0 0 ];
% maxr=450;
% plot_tuning2d(eod(2),file(2),18,'t_col',6,'clim',[-2 2],...
%     'x_lim',[-maxr maxr],'y_lim',[-maxr maxr],...
%     'maxr',maxr,'minsamp',5,...
%     'image',IMGS.cdata,'mfunc','offset','ops',ops,'bgcol',bgc,'upsamp',1.5);
% set(gca,'Alim',[1 12]);

get_tail_wall(ha([4 7]));
set(ha(4),'Xlim',[-28 28],'Ylim',[-23 23],'Clim',[-0.05 0.03]);
axis(ha(4),'image');
set(ha(7),'Xlim',[-28 28],'Ylim',[-23 23],'Clim',[-0.025 0.025]);
axis(ha(7),'image');

axes(ha(8));
ha(8).clo;
[Mz,N0,al]=plot_tuning_polar(eod(2),file(2),18,'t_col',6,'clim',[-1 1],...
    'image',str_image.cdata,'mfunc','slope','ops',ops,'bgcol',bgc,...
    'upsamp',0.25,'r_nbins',15,'th_nbins',15,'poly',p);
ha(8).Position=pos{8}-[0.055 0 0 0];
% plot_tuning2d(eod(2),file(2),18,'t_col',6,'clim',[-1 1],...
%     'x_lim',[-maxr maxr],'y_lim',[-maxr maxr],...
%     'maxr',maxr,'minsamp',1,'x_nbins',20,'y_nbins',20,...
%     'image',IMGS.cdata,'mfunc','slope','ops',ops,'bgcol',bgc,'upsamp',1.5);
ha(3).Position=pos{3}+[0.05 0 -0.05 0];
[hlfp,plfp]=split_axes(ha(3),3,2,0.01,0.01);
for i=1:numel(al)
    hlfp(i).clo;
%     set(ca(i),'XLim',[-10 40]);    
    copyaxes(al(i),hlfp(i),1);
    set(hlfp(i),'Colormap',al(i).Colormap,'Color','none');
end

%% model figures
% ind=find(r>20 & r<40 & ~isnan(t));
% [mm,mind]=min(th(ind));
% [fidx,fridx]=get_frame_indices(ind(mind),frame(2).ind0);
% filenum=num2str(ops.seg(2).files(fidx));
% vidname = [getenv('DATAPATH'),'\',getenv('SESSDATE'),'\','video_',filenum];
% vid=VideoReader([vidname,'.avi']);    
% trackfile=[vidname,ops.trackname];
% [num,txt,raw] = xlsread([trackfile,'.csv']);    
% [FF,ax_grid]=visualize_fish_modeling(vid,txt,num,fridx{1},1,file(2).circle);
% 
% ind=find(r>180 & ~isnan(t));
% [mm,mind]=min(t(ind));
% [MM,Mind]=max(t(ind));
% [Mm,Mmind]=min(abs(t(ind)-median(t)));
% [fidx,fridx]=get_frame_indices(ind(mind),frame(2).ind0);
% filenum=num2str(ops.seg(2).files(fidx));
% vidname = [getenv('DATAPATH'),'\',getenv('SESSDATE'),'\','video_',filenum];
% vid=VideoReader([vidname,'.avi']);    
% trackfile=[vidname,ops.trackname];
% [num,txt,raw] = xlsread([trackfile,'.csv']);    
% [FF,ax_grid]=visualize_fish_modeling(vid,txt,num,fridx{1},1,file(2).circle,file(2).BG);
% 
% [fidx,fridx]=get_frame_indices(ind(Mind),frame(2).ind0);
% filenum=num2str(ops.seg(2).files(fidx));
% vidname = [getenv('DATAPATH'),'\',getenv('SESSDATE'),'\','video_',filenum];
% vid=VideoReader([vidname,'.avi']);    
% trackfile=[vidname,ops.trackname];
% [num,txt,raw] = xlsread([trackfile,'.csv']);    
% [FF,ax_grid]=visualize_fish_modeling(vid,txt,num,fridx{1},1,file(2).circle,file(2).BG);
% 
% [fidx,fridx]=get_frame_indices(ind(Mmind),frame(2).ind0);
% filenum=num2str(ops.seg(2).files(fidx));
% vidname = [getenv('DATAPATH'),'\',getenv('SESSDATE'),'\','video_',filenum];
% vid=VideoReader([vidname,'.avi']);    
% trackfile=[vidname,ops.trackname];
% [num,txt,raw] = xlsread([trackfile,'.csv']);    
% [FF,ax_grid]=visualize_fish_modeling(vid,txt,num,fridx{1},1,file(2).circle,file(2).BG);
%% ex-affernce / re-afference stats
load('Z:\analysis_data\ELL__Models_Database.mat');
axes(ha(6));
ha(6).clo;
[ rsq]= plot_tuning_exafference(eod(2),file(2),18,'t_col',6,'ind_lim',1e5+[1 1e4])
set(gca,'Xlim',[-2 2],'Ylim',[-1 1],'FontSize',16,'YTick',[-1 -.5 0 .5 1],'YTickLabel',{'-1','-.5','0','.5','1'})
%inset
hd=axes();
hd.Position=pos{6}.*[1 1 .1 .4]+[pos{6}(3)/10 pos{6}(4)*.55 0 0];
hd.clo;
boxplot(groups.EX_Model,ones(size(groups.EX_Model)),'PlotStyle','compact','Colors',0.5*[1 1 1]);
bx=findobj('Tag','Box');
set(bx,'LineWidth',15);
ws=findobj('Tag','Whisker');
set(ws,'LineWidth',3);
cro=findobj('Tag','MedianInner');
set(cro,'Marker','none');
cri=findobj('Tag','MedianOuter');
set(gca,'YLim',[0 .35],'YTick',[0 .1 .2 .3],'YTickLabel',{'0','.1','.2','.3'},'XTickLabel',[],'FontSize',12);
% R=[];
% for i=1:numel(database.experiments)
%     for j=1:numel(database.experiments(i).groups);
%         R=[R;database.experiments(i).groups(j).model(1).R];
%     end
% end
% violin(R,'facecolor',0.5*[1 1 1],'edgecolor',0.25*[1 1 1],'facealpha',1,'mc','','medc','');
% set(gca,'FontSize',12,'XTick',[],'Ylim',[0 .4]);

% R=struct2table(database.r2);
% M=[mean(R.ex) mean(R.re) mean(R.all) ];
% Mn=[mean(R.all-R.ex_ctl) mean(R.all-R.re_ctl)];
% S=[std(R.ex) std(R.re) std(R.all)];
% Sn=[std(R.all-R.ex_ctl) std(R.all-R.re_ctl)];
% figure;
% B=bar([1 2 3],M);
% set(B,'FaceColor',0.5*[1 1 1],'LineStyle','none');
% hold on;
% E=errorbar([1 2 3],M,S);
% set(E,'LineStyle','none','LineWidth',3,'Color',[0 0 0]);
% Bn=bar([2 3],-Mn);
% set(Bn,'FaceColor',1*[1 1 1],'LineStyle','-');
% En=errorbar([2 3],-Mn,Sn);
% set(En,'LineStyle','none','LineWidth',3,'Color',[0 0 0]);
% set(gcf,'Position', [464    88   200   195]);

%% reafference variance
% load('Z:\analysis_data\database_LFP_wall.mat');
COL=brewermap(12,'Paired');
% reafference non uniformity
axes(ha(9));
S=scatter(groups.Tail_Var,groups.Tail_Var_Ctl,100,COL(2,:),'filled');
set(S,'MarkerFaceAlpha',0.5);
hold on;
S=scatter(groups.Rate_Var,groups.Rate_Var_Ctl,100,COL(8,:),'filled','square');
set(S,'MarkerFaceAlpha',0.5);
H=plot([0 0.2],[0 0.2],':k');
set(gca,'FontSize',16,'Xlim',[0 0.17],'Ylim',[0 0.17],'YTick',[0 .05 .1 .15],'YTickLabel',{'0','.05','.1','.15'},'Color','none');

%inset
hh=axes;
hh.Position=pos{9}.*[1 1 .4 .4]+[pos{9}(3)/6 pos{9}(4)*.55 0 0];

X=[groups.Tail_Var';groups.Tail_Var_Ctl';groups.Rate_Var';groups.Rate_Var_Ctl']';
% G=reshape(ones(size(groups,1),1)*[1 2 3 4],[],1);
% boxplot(X,G,'PlotStyle','compact');
% violin(X,'facecolor',COL([2 1 8 7],:),'edgecolor',1*[1 1 1],'facealpha',1,'mc','','medc','');
% legend('off');
M=mean(X);
S=std(X);
CC=COL([2 1 8 7],:);
for i=1:size(X,2)
    B=bar(i,M(i));
    set(B,'FaceColor',CC(i,:),'LineStyle','none');
    hold on;
    E=errorbar(i,M(i),S(i));
    set(E,'LineStyle','none','LineWidth',3,'Color',CC(i,:)*.5);
end
set(hh,'FontSize',12,'XTick',[],'Ylim',[0 0.2],'YTick',[0 0.2],'YAxisLocation','left')
set(gcf,'Position', [148    0  1200   1000],'Color','none');

