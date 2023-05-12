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
load('ROIs.mat','p_IO_new');
%% plot ex afference, re afference maps
Cslope=[-1 1];
bgcol=[1 1 1];
upsamp=1.5;
minsamp=10;
nline=0.25;
wline=0.5;
fontsize=7;
y_nbins=30;
COL=[236 157 118;210 149 0;74 133 34]/255;
unit1_cmap=ones(4,1)*COL(2,:);
unit2_cmap=ones(4,1)*COL(3,:);

unit1_eod=eod2(6); unit1_file=file2(6); unit1_col=24; unit1_ind=ops2.seg(6).ind_lim(2,:); unit1_clim=[0 80]; unit1_lfpcol=21; unit1_rast=unit1_eod.raster{2};
unit2_eod=eod1(4); unit2_file=file1(4); unit2_col=20; unit2_ind=ops1.seg(4).ind_lim(2,:); unit2_clim=[0 80]; unit2_lfpcol=17; unit2_rast=unit2_eod.raster{2};

Fall=figure;
[ha, pos] = tight_subplot(3, 2, [.0 .02],[.0 0.01],[.1 0.0],[0.45 0.275 0.275],[]); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
set(Fall,'Units','centimeters');
set(Fall,'Position', [2  2  12  12],'Color',[1 1 1]);
ha(1).Position=pos{1}+[0.05 0.05 -0.1 -0.05];
ha(2).Position=pos{2}+[0.05 0.05 -0.1 -0.05];
axes(ha(3));
[M1,N1,s1]=plot_tuning2d(unit1_eod,unit1_file,unit1_col,'mfunc','offset','t_col',unit1_lfpcol,'image',str_image.cdata,...
    'clim',unit1_clim,'bgcol',bgcol,'x_lim',[-550 550],'y_lim',[-500 300],'y_nbins',y_nbins,...
    'minsamp',minsamp,'upsamp',upsamp,'zscale',[0 40],'ind_lim',unit1_ind,...
    'rastind',2,'ops',ops1,'roinum',numel(p_IO),'poly',p_IO_new,'aux_axes',ha(1),'cmaps',unit1_cmap,'fit_mode','intercept','psth_lines',3);
axes(ha(5));
plot_tuning2d(unit1_eod,unit1_file,unit1_col,'mfunc','slope','t_col',unit1_lfpcol,'image',str_image.cdata,...
    'clim',Cslope,'bgcol',bgcol,'x_lim',[-550 550],'y_lim',[-500 300],'y_nbins',y_nbins,...
    'minsamp',minsamp,'upsamp',upsamp,'ind_lim',unit1_ind);

axes(ha(4));
plot_tuning2d(unit2_eod,unit2_file,unit2_col,'mfunc','offset','t_col',unit2_lfpcol,'image',str_image.cdata,...
    'clim',unit2_clim,'bgcol',bgcol,'x_lim',[-550 550],'y_lim',[-550 300],'y_nbins',y_nbins,...
    'minsamp',minsamp,'upsamp',upsamp,'zscale',[0 40],'ind_lim',unit2_ind,...
    'rastind',2,'ops',ops2,'roinum',numel(p_IO),'poly',p_IO_new,'aux_axes',ha(2),'cmaps',unit2_cmap,'fit_mode','intercept','psth_lines',3);
axes(ha(4));
colorbar('off');
axes(ha(6));
plot_tuning2d(unit2_eod,unit2_file,unit2_col,'mfunc','slope','t_col',1,'image',str_image.cdata,...
    'clim',Cslope,'bgcol',bgcol,'x_lim',[-550 550],'y_lim',[-550 300],'y_nbins',y_nbins,...
    'minsamp',minsamp,'upsamp',upsamp,'ind_lim',unit2_ind);
colorbar('off');

for i=1:numel(p_IO_new)
    P(i,:)=[min(p_IO_new{i}) range(p_IO_new{i})];
    rectangle(ha(3),'Position',P(i,:),'Curvature',[1 1]);
    rectangle(ha(4),'Position',P(i,:),'Curvature',[1 1]);
    rectangle(ha(5),'Position',P(i,:),'Curvature',[1 1]);
    rectangle(ha(6),'Position',P(i,:),'Curvature',[1 1]);
end
C=findobj('Type','Colorbar');
set(C(:),'Location','west','AxisLocation','out','LineWidth',nline,'FontSize',fontsize);
for i=1:numel(C)
    C(i).Position(3)=C(i).Position(3)-0.02;
    C(i).Position(1)=C(i).Position(1)-0.02;
end
%% export
set(Fall,'Color','none');
hall=findobj();
for i=2:numel(hall)
    hall(i).Visible='on';
end

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


