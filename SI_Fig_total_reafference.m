% make panels for Figure 1: Afferent input in freely moving fish
[ret,sname]=system('hostname');
if(strfind(sname,'Mac'))
    rootdir='/Volumes/aw3057/';
    sep='/';
else
    rootdir='Z:\';
    sep='\';
end
COL=[226 139 138;210 149 0;74 133 34]/255;

CHSL=rgb2hsl(COL);
nline=0.25;
wline=0.5;
bwidth=4;
msize=.5;
fontsize=7;
%% load data
% path(path,'C:\Users\Ephys\iCloudDrive\Documents\Mormyrid_Data\electric model\new');
% path(path,'Z:\GitHub\Fish-Pipeline')
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

%% plot ex afference, re afference maps
%inds=[5e5 inf];
Coff=[-1.5 1.5];
Cslope=[-.5 .5];
% bgcol=[0 0 0];
bgcol=[1 1 1];
upsamp=.25;
minsamp=20;

unit1_cmap=ones(3,1)*COL(2,:);
unit2_cmap=ones(3,1)*COL(3,:);

r_nbins=10;
th_nbins=15;
% ylim=[-550 250];
% xlim=[-550 550];
% unit1_eod=eod1(4); unit1_file=file1(4); unit1_col=19; unit1_ind=ops1.seg(4).ind_lim(1,:); unit1_clim=[0 2]; unit1_lfpcol=17; 
% unit1_ops=ops1; unit1_rastind=1; unit1_rast=unit1_eod.raster{1};
% unit2_eod=eod1(4); unit2_file=file1(4); unit2_col=20; unit2_ind=ops1.seg(4).ind_lim(2,:); unit2_clim=[0 2]; unit2_lfpcol=17; 
% unit2_ops=ops1; unit2_rastind=2; unit2_rast=unit2_eod.raster{2};
unit1_eod=eod2(6); unit1_file=file2(6); unit1_col=24; unit1_ind=ops2.seg(6).ind_lim(2,:); unit1_clim=[0 2]; unit1_lfpcol=21; 
unit1_ops=ops2; unit1_rastind=2; unit1_rast=unit1_eod.raster{2};lfp_chan=13;
unit2_eod=eod2(6); unit2_file=file2(6); unit2_col=23; unit2_ind=ops2.seg(6).ind_lim(1,:); unit2_clim=[0 1.5]; unit2_lfpcol=21; 
unit2_ops=ops2;  unit2_rastind=1; unit2_rast=unit2_eod.raster{1};
% lfp_eod=eod1(2); lfp_file=file1(2); lfp_col=18; lfp_ops=ops1; 
%% get ROIs 
rois=0;
if(rois)
    alim=[-1 3];
    r_nbins=10;
    th_nbins=15;
    upsamp=1;
    
    figure;
    h1=subplot(2,2,1); h2=subplot(2,2,2); 
    h3=subplot(2,2,3); h4=subplot(2,2,4); 

    axes(h3);
    [M2,N2,hh]=plot_tuning_polar(unit2_eod,unit2_file,unit2_col,'mfunc','slope','t_col',unit2_lfpcol,'image',str_image.cdata,...
        'clim',Cslope,'bgcol',bgcol,'r_nbins',r_nbins,'th_nbins',th_nbins,'alim',alim,...
        'minsamp',minsamp,'upsamp',upsamp,'ind_lim',unit2_ind,'cmaps',unit2_cmap);
%         'minsamp',minsamp,'upsamp',upsamp,'ind_lim',unit2_ind,'rastind',unit2_rastind,'roinum',3,'ops',unit2_ops,'aux_axes',h4,'cmaps',unit2_cmap);

    axes(h1);    
    [M1,N1,hh]=plot_tuning_polar(unit1_eod,unit1_file,unit1_col,'mfunc','slope','t_col',unit1_lfpcol,'image',str_image.cdata,...
        'clim',Cslope,'bgcol',bgcol,'r_nbins',r_nbins,'th_nbins',th_nbins,'alim',alim,...
        'minsamp',minsamp,'upsamp',upsamp,'ind_lim',unit1_ind,'rastind',unit1_rastind,'roinum',3,'ops',unit1_ops,'aux_axes',h2,'cmaps',unit1_cmap);
%         'minsamp',minsamp,'upsamp',upsamp,'ind_lim',unit1_ind,'cmaps',unit1_cmap);
    p_SI=hh.p;    
    axes(h3);
    [M2,N2,hh]=plot_tuning_polar(unit2_eod,unit2_file,unit2_col,'mfunc','slope','t_col',unit2_lfpcol,'image',str_image.cdata,...
        'clim',Cslope,'bgcol',bgcol,'r_nbins',r_nbins,'th_nbins',th_nbins,'alim',alim,...
        'minsamp',minsamp,'upsamp',upsamp,'ind_lim',unit2_ind,'rastind',unit2_rastind,'poly',p_SI,'ops',unit2_ops,'aux_axes',h4,'cmaps',unit2_cmap);

    save('scatter_polygons.mat','p_SI','-append');
%     axes(h1);
%     [M1,N1,hh]=plot_tuning_polar(unit1_eod,unit1_file,unit1_col,'mfunc','slope','t_col',6,'image',str_image.cdata,...
%         'clim',Cslope,'bgcol',bgcol,'r_nbins',r_nbins,'th_nbins',th_nbins,'alim',alim,...
%         'minsamp',minsamp,'upsamp',upsamp,'ind_lim',unit1_ind,'rastind',unit1_rastind,'poly',p,'ops',unit1_ops,'aux_axes',h2,'cmaps',unit1_cmap);
    
else
    load('scatter_polygons.mat');
end

%%
% alim=[-1 3];
r_nbins=10;
th_nbins=10;
upsamp=.25;
binsize=1;

Fall=figure;
set(Fall,'Units','centimeters');
set(Fall,'OuterPosition', [2  2  11.4  15],'Color',bgcol);
[h0, pos0] = tight_subplot(3, 1, [.075 .0],[.035 0.02],[.05 0.0],[.35 .35 .3]); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
%% first row - type 1 
[ha,pa]=tight_subplot(1, 3, [.0 .02],[.0 .0],[.00 .05],[],[.5 .125 .375],h0(1)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]

axes(ha(3));
ha(3).clo;
[M1,N1,hh]=plot_tuning_polar(unit1_eod,unit1_file,unit1_col,'mfunc','slope','t_col',unit1_lfpcol,'image',str_image.cdata,...
    'clim',Cslope,'bgcol',bgcol,'r_nbins',r_nbins,'th_nbins',th_nbins,'binsize',binsize,...
    'minsamp',minsamp,'upsamp',upsamp,'ind_lim',unit1_ind,'rastind',unit1_rastind,'poly',p_SI,'ops',unit1_ops,'aux_axes',ha(1),'cmaps',unit1_cmap);

axes(ha(2));
m=M1(:);
n=N1(:);
n=n/max(n)*50+1;
S=scatter(randn(size(m))*.1,m,n,m,'filled','MarkerEdgeColor',0.25*[1 1 1],'MarkerFaceAlpha',1,'LineWidth',.1);
set(ha(2),'Xlim',[-1 1],'Ylim',Cslope,'YTickLabel',[],'Clim',Cslope,'XAxisLocation','origin','XTick',[]);
set(ha(2),'Colormap',ha(3).Colormap,'Color','none');

%% second row - type 2 
[hb,pb]=tight_subplot(1, 3, [.0 .02],[.0 .0],[.00 .05],[],[.5 .125 .375],h0(2)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]

axes(hb(3));
hb(3).clo;
[M2,N2,hh]=plot_tuning_polar(unit2_eod,unit2_file,unit2_col,'mfunc','slope','t_col',unit2_lfpcol,'image',str_image.cdata,...
    'clim',Cslope,'bgcol',bgcol,'r_nbins',r_nbins,'th_nbins',th_nbins,'binsize',binsize,...
    'minsamp',minsamp,'upsamp',upsamp,'ind_lim',unit2_ind,'rastind',unit2_rastind,'poly',p_SI,'ops',unit2_ops,'aux_axes',hb(1),'cmaps',unit2_cmap);

axes(hb(2));
m=M2(N2>minsamp);
n=N2(N2>minsamp);
n=n/max(n)*50+1;
S=scatter(randn(size(m))*.1,m,n,m,'filled','MarkerEdgeColor',0.25*[1 1 1],'MarkerFaceAlpha',1,'LineWidth',.1);
set(hb(2),'Xlim',[-1 1],'Ylim',Cslope,'YTickLabel',[],'Clim',Cslope,'XAxisLocation','origin','XTick',[]);
set(hb(2),'Colormap',hb(3).Colormap,'Color','none');

C=findobj('Type','Colorbar');
set(C(:),'Location','east','AxisLocation','out','LineWidth',nline,'TickLength',0.025,'TickDirection','both','FontSize',fontsize,'FontWeight','normal');
for i=1:numel(C)
    C(i).Position(3)=C(i).Position(3)*.75;
    C(i).Position(1)=C(i).Position(1)+0.05;
end

%% third row
h0(3).Visible='off';

%% expand figure, remove all panels save for images for bitmap copying
set(Fall,'Color','none');
ha(3).Color='none';
hb(3).Color='none';
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

