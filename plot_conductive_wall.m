% make panels for Figure 1: Afferent input in freely moving fish
%% load data
path(path,'Z:\GitHub\Fish-Model');
path(path,'Z:\GitHub\Fish-Pipeline')
if(~exist('eod'))
    load('Z:\analysis_data\20221116\output.mat')
end
load('Z:\mormyrid_data\fish_images\fish_silhouette.mat');     
str_image=IMGS;
setenv('PXLSIZE','0.423'); %pxl to mm conversion
pxl2mm=0.423;
setenv('FRAMESCALE',num2str(ops.framescale));
setenv('DATAPATH','Z:\mormyrid_data');
%% start figure
Fall=figure;
bgc=[0 0 0];
% bgc=[1 1 1];
fontsize=16;
[ha, pos] = tight_subplot(2, 4, [.02 .01],[.0 .0],[.0 .0]);
% ha(1).Visible='off'; %place for moodel scheme

%% plot ex afference, re afference maps
get_tail_wall(ha([1 2]),0.65*pi,1);

axes(ha(3));
ha(3).clo;
[Ex1,N01]=plot_tuning_polar(eod(1),file(1),[15:18],'zfunc',@nanmean,'t_col',6,'clim',[-2 2],...
    'image',IMGS.cdata,'mfunc','offset','ops',ops,'bgcol',bgc,...
    'upsamp',0.25,'r_nbins',10,'th_nbins',15);
set(ha(3),'CLim',[-1 1.5]);
colorbar('off');

axes(ha(4));
ha(4).clo;
[Re1,N01]=plot_tuning_polar(eod(1),file(1),[15:18],'zfunc',@nanmean,'t_col',6,'minsamp',25,...
    'image',IMGS.cdata,'mfunc','slope','ops',ops,'bgcol',bgc,...
    'upsamp',0.25,'r_nbins',10,'th_nbins',15);
set(ha(4),'CLim',[-.35 .35]);
colorbar('off');

get_tail_wall(ha([5 6]),0.65*pi,-1);

axes(ha(7));
ha(7).clo;
[Ex2,N02]=plot_tuning_polar(eod(2),file(2),[15:18],'zfunc',@nanmean,'t_col',6,'minsamp',25,...
    'image',IMGS.cdata,'mfunc','offset','ops',ops,'bgcol',bgc,...
    'upsamp',0.25,'r_nbins',10,'th_nbins',15);
set(ha(7),'CLim',[-1 1]);
colorbar('off');

axes(ha(8));
ha(8).clo;
[Re2,N02]=plot_tuning_polar(eod(2),file(2),[15:18],'zfunc',@nanmean,'t_col',6,'minsamp',25,...
    'image',IMGS.cdata,'mfunc','slope','ops',ops,'bgcol',bgc,...
    'upsamp',0.25,'r_nbins',10,'th_nbins',15);
set(ha(8),'CLim',[-.35 .35]);
colorbar('off');


