% make video showing tail effects on response to brass pole for example
% LFP, type 1 and type 2 cells
[ret,sname]=system('hostname');
if(strfind(sname,'Mac'))
    rootdir='/Volumes/aw3057/';
    sep='/';
else
    rootdir='Z:\';
    sep='\';
end
COL=[226 139 138;210 149 0;74 133 34]/255;
nline=0.25;
wline=0.5;
bwidth=4;
msize=.5;
fontsize=7;
%% load data
if(~exist('eod1'))
    load([rootdir,'analysis_data',sep,'20190131',sep,'output.mat']);
    eod1=eod;
    file1=file;
    ops1=ops;
end
if(~exist('eod2'))
    load([rootdir,'analysis_data',sep,'20190624',sep,'output.mat']);
    eod2=eod;
    file2=file;
    ops2=ops;
end

load('Z:\mormyrid_data\fish_images\fish_silhouette.mat');
str_image=IMGS;
setenv('PXLSIZE','0.423'); %pxl to mm conversion
pxl2mm=0.423;
setenv('FRAMESCALE',num2str(ops.framescale));
setenv('SESSDATE','20190131');
setenv('DATAPATH','Z:\mormyrid_data');

seg=6; fnum=9; ss=13800; 
fdate=num2str(ops2.seg(seg).dates(1));
dname=[ops.datapath,'\',fdate,'\'];
fname=ops.seg(seg).files(fnum); %file to take
framenum=200; %number of frames
frameind=ss+[0:framenum];
if(~exist('fish_frames'))
    vidname=[dname,'video_',num2str(fname)];
    trackfile=[dname,'video_',num2str(fname),ops.trackname];
    vid = VideoReader([vidname,'.avi']);
    F=figure;
    set(F,'Color',[0 0 0],'Position',[1 41 1920 964]);
    setenv('DATAPATH',ops.datapath);
    setenv('SESSDATE',fdate);
    setenv('FRAMESCALE',num2str(ops.framescale));
    [fish_frames,ax_grid] = get_fish_image(fname,frameind,'headfix','off','track','off');
%     [frames,ax_grid]=visualize_fish_tracking(vidname,trackfile,ss,framenum)
end

%% generate fish images
K=51;
r = @(x) sqrt(x(:,1).^2 + x(:,2).^2);
w = @(x) atan2(x(:,2), x(:,1));
I0=str_image.cdata;
a=linspace(-1.35,1.35,K);
figure;
for i=1:numel(a)
    q = @(x) w(x) - a(i)*(x(:,2)>.0) .* (r(x).^2);
    f = @(x) [(r(x)) .* cos(q(x)), (r(x)) .* sin(q(x))];
    g = @(x, unused) f(x);

    tform2 = maketform('custom', 2, 2, [], g, []);
    I{i} = imtransform(I0, tform2, 'UData', [-1 1], 'VData', [-1 1]+.1, ...
        'XData', [-1 1], 'YData', [-1 1],'FillValues',255);
    imshow(I{i});
    drawnow;
end
%% 
Coff=[-1 1];
Cslope=[-.5 .5];
bgcol=[0 0 0];
upsamp=2;
minsamp=3;
lfp_eod=eod1(4); lfp_file=file1(4); lfp_col=17; %lfp_zscale=[430 40];
unit1_eod=eod2(6); unit1_file=file2(6); unit1_col=24; unit1_ind=ops2.seg(6).ind_lim(2,:); unit1_clim=[0 2]; unit1_lfpcol=21; unit1_rast=unit1_eod.raster{2};
unit2_eod=eod1(4); unit2_file=file1(4); unit2_col=20; unit2_ind=ops1.seg(4).ind_lim(2,:); unit2_clim=[0 2]; unit2_lfpcol=17; unit2_rast=unit2_eod.raster{2};

x_nbins=23;
y_nbins=23;

p=linspace(0,1,K+1);
% p=p(1:end-1);
p0=.35;
p1=linspace(0,1-p0,K);
p2=linspace(p0,1,K);
p3=[p1;p2];


% q=quantile(unit1_eod.data(:,6),p);
%% LFP
% q=reshape(quantile(lfp_eod.data(:,6),p3(:)),2,[]);
t=lfp_eod.data(:,6);
r0=range(t)*p0;
q1=linspace(min(t),max(t)-r0,K);
q2=linspace(min(t)+r0,max(t),K);
q3=[q1;q2];

figure;
for i=1:K
    plot_tuning2d(lfp_eod,lfp_file,lfp_col,'mfunc','mean','t_col',6,'t_val',q3(:,i)',...
        'image',I{i},'clim',[-1 1],'bgcol',bgcol,'x_lim',[-550 550],'y_lim',[-550 250],...
        'x_nbins',x_nbins,'y_nbins',y_nbins,'zscore',0,'minsamp',minsamp,'upsamp',upsamp);
    colorbar('off');
    set(gca,'Alim',[0 12]);
    frame_lfp(i)=getframe;
end

plot_tuning2d(lfp_eod,lfp_file,lfp_col,'mfunc','mean','t_col',6,'image',str_image.cdata,...
    'clim',Coff,'bgcol',bgcol,'x_lim',[-550 550],'y_lim',[-550 250],...
    'minsamp',20,'upsamp',upsamp);
colorbar('off');
frame_lfp_offset=getframe;

plot_tuning2d(lfp_eod,lfp_file,lfp_col,'mfunc','slope','t_col',6,'image',str_image.cdata,...
    'clim',Cslope,'bgcol',bgcol,'x_lim',[-550 550],'y_lim',[-550 250],...
    'minsamp',20,'upsamp',upsamp);
colorbar('off');
frame_lfp_slope=getframe;
%% unit 1
t=unit1_eod.data(:,6);
r0=range(t)*p0;
q1=linspace(min(t),max(t)-r0,K);
q2=linspace(min(t)+r0,max(t),K);
q3=[q1;q2];

% q=reshape(quantile(unit1_eod.data(:,6),p3(:)),2,[]);
figure;
for i=1:K
    plot_tuning2d(unit1_eod,unit1_file,unit1_col,'mfunc','mean','t_col',6,'t_val',q3(:,i)',...
        'image',I{i},'clim',[0 2],'bgcol',bgcol,'x_lim',[-550 550],'y_lim',[-550 250],...
        'ind_lim',unit1_ind,...
        'x_nbins',x_nbins,'y_nbins',y_nbins,'zscore',0,'minsamp',minsamp,'upsamp',upsamp);
    colorbar('off');
    set(gca,'Alim',[0 12]);
    frame_unit1(i)=getframe;
end

plot_tuning2d(unit1_eod,unit1_file,unit1_col,'mfunc','mean','t_col',6,'image',str_image.cdata,...
    'clim',unit1_clim,'bgcol',bgcol,'x_lim',[-550 550],'y_lim',[-500 300],...
    'minsamp',20,'upsamp',upsamp,'zscore',0,'ind_lim',unit1_ind);
colorbar('off');
frame_unit1_offset=getframe;

plot_tuning2d(unit1_eod,unit1_file,unit1_col,'mfunc','slope','t_col',6,'image',str_image.cdata,...
    'clim',Cslope,'bgcol',bgcol,'x_lim',[-550 550],'y_lim',[-500 300],...
    'minsamp',20,'upsamp',upsamp,'ind_lim',unit1_ind);
colorbar('off');
frame_unit1_slope=getframe;

%% unit 2
t=unit2_eod.data(:,6);
r0=range(t)*p0;
q1=linspace(min(t),max(t)-r0,K);
q2=linspace(min(t)+r0,max(t),K);
q3=[q1;q2];

% q=reshape(quantile(unit2_eod.data(:,6),p3(:)),2,[]);
figure;
for i=1:K
    plot_tuning2d(unit2_eod,unit2_file,unit2_col,'mfunc','mean','t_col',6,'t_val',q3(:,i)',...
        'image',I{i},'clim',[0 2],'bgcol',bgcol,'x_lim',[-550 550],'y_lim',[-550 250],...
        'ind_lim',unit2_ind,...
        'x_nbins',x_nbins,'y_nbins',y_nbins,'zscore',0,'minsamp',minsamp,'upsamp',upsamp);
    colorbar('off');
    set(gca,'Alim',[0 12]);
    frame_unit2(i)=getframe;
end

plot_tuning2d(unit2_eod,unit2_file,unit2_col,'mfunc','mean','t_col',6,'image',str_image.cdata,...
    'clim',unit2_clim,'bgcol',bgcol,'x_lim',[-550 550],'y_lim',[-550 250],...
    'minsamp',20,'upsamp',upsamp,'zscore',0,'ind_lim',unit2_ind);
colorbar('off');
frame_unit2_offset=getframe;

plot_tuning2d(unit2_eod,unit2_file,unit2_col,'mfunc','slope','t_col',6,'image',str_image.cdata,...
   'clim',Cslope,'bgcol',bgcol,'x_lim',[-550 550],'y_lim',[-550 250],...
    'minsamp',20,'upsamp',upsamp,'ind_lim',unit2_ind);
colorbar('off');
frame_unit2_slope=getframe;

%% compose video
frate=round(K/2);
d_2s=linspace(0,1,2*frate);
d_1s=linspace(0,1,frate);
d_05s=linspace(0,1,frate/2);

vidpath='C:\Users\sawte\Dropbox\efish_paper\Rev1\Videos\';
vidout=VideoWriter([vidpath,'WallachSawtell2023_MovieS3'],'mpeg-4');
vidout.FrameRate=frate;
vidout.Quality=10;
open(vidout);
Fall=figure;
set(Fall,'units','pixels','Position',[100 100 492*2 276*2],'Color',bgcol)
%% title page
title_image=imread('C:\Users\sawte\Dropbox\efish_paper\Videos\pole_vid_title.png');  
A=axes('Position',[0 .2 1 .8],'Color',bgcol);
imshow(title_image);
A2=axes('Position',[0.1,0.1,0.43,0.47],'Color',bgcol);
Ta=annotation('textarrow','String','Brass Object','FontSize',14,'X',[0.28 0.31],'Y',[0.4 0.33],'Color',[1 1 1]);
Abox=axes('Position',[0 0 1 1],'Color','none','XColor','none','YColor','none');
R=rectangle('Position',[0 0 1 1],'LineStyle','none','FaceColor',[0 0 0 1]);
L0=annotation('textbox','String','0.5X Speed','Color',[1 1 1],'FontSize',14,'Position',[0.1 0.01 0.2 0.1],'LineStyle','none');
% k=1;
M=frate*4;
for i=1:M
    axes(A2);
    Im=fish_frames(i).cdata;
    x0=413; y0=415; Rad=400;
    [X,Y]=meshgrid(1:size(Im,2),1:size(Im,1));
    blnkind=find(((X-x0).^2 + (Y-y0).^2)>Rad^2);
    CD=flipud(Im(:,:,1));
    CD(blnkind)=0;
    Im=repmat(CD,1,1,3);    
    imshow(Im);
    D=numel(d_05s);
    k=max([D-i+1,1,i-(M-D)]);
    set(Ta,'Color',(1-d_05s(k))*[1 1 1]);        
    axes(Abox);
    R.FaceColor=[0 0 0 d_05s(k)];
    writeVideo(vidout,getframe(Fall));    
end
delete(A);
delete(A2);
delete(Ta);
delete(Abox);
delete(L0);
%%
[ha, pos] = tight_subplot(2,3, [.0 .0],[.15 0.15],[.1 .1]); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
set(ha(:),'Color',bgcol,'Xcolor','none','Ycolor','none');
lfp_text='LFP (electrosensory input)';
unit1_text='Type 1 (putative interneuron)';
unit2_text='Type 2 (putative output)';
%% LFP movie
Abox=axes('Position',pos{1}.*[1 1 1 1.2],'Color','none','XColor','none','YColor','none');
R=rectangle('Position',[0 0 1 1],'LineStyle','none','FaceColor',[0 0 0 1]);
k=1;
for j=1:2
    for i=[1:(2*K-1)]+floor(K/2)
        axes(ha(1));
        frk=abs(K-1-abs(i-K))+1;
        imshow(frame_lfp(frk).cdata);
        T=title(lfp_text,'FontSize',14,'FontName','Arial','FontWeight','bold','Color',COL(1,:));        
        
        if(k<=D)
            axes(Abox);    
            R.FaceColor=[0 0 0 1-d_05s(k)];
            k=k+1;
        end
        writeVideo(vidout,getframe(Fall));
    end
end

% LFP maps
Abox=axes('Position',pos{4}.*[0 1 2 1],'Color','none','XColor','none','YColor','none');
R=rectangle('Position',[0 0 1 1],'LineStyle','none','FaceColor',[0 0 0 1]);

axes(ha(1));
ha(1).clo;
imshow(frame_lfp_offset.cdata);
T=title(lfp_text,'FontSize',14,'FontName','Arial','FontWeight','bold','Color',COL(1,:));        
ylabel(ha(1),{'Average response','(tuning map'},'Color',[1 1 1]...
    ,'FontSize',12,'FontName','Arial','FontWeight','bold','Visible','on');
axes(ha(4));
imshow(frame_lfp_slope.cdata);
ylabel(ha(4),{'Tail regression map','(object coordinates)'},'Color',[1 1 1]...
    ,'FontSize',12,'FontName','Arial','FontWeight','bold','Visible','on');

for i=1:2*frate   
    if(i<=D)
        axes(Abox);        
        R.FaceColor=[0 0 0 1-d_05s(i)];
    end
    writeVideo(vidout,getframe(Fall));
end
%% unit 1 movie
Abox=axes('Position',pos{2}.*[1 1 1 1.2],'Color','none','XColor','none','YColor','none');
R=rectangle('Position',[0 0 1 1],'LineStyle','none','FaceColor',[0 0 0 1]);
k=1;

for j=1:2
    for i=[1:(2*K-1)]+floor(K/2)
        axes(ha(2));
        frk=abs(K-1-abs(i-K))+1;        
        imshow(frame_unit1(frk).cdata);
        T=title(unit1_text,'FontSize',14,'FontName','Arial','FontWeight','bold','Color',COL(2,:));
        if(k<=D)
            axes(Abox);    
            R.FaceColor=[0 0 0 1-d_05s(k)];
            k=k+1;
        end        
        writeVideo(vidout,getframe(Fall));
    end
end

Abox=axes('Position',pos{5}.*[1 1 1 1],'Color','none','XColor','none','YColor','none');
R=rectangle('Position',[0 0 1 1],'LineStyle','none','FaceColor',[0 0 0 1]);
axes(ha(2));
ha(2).clo;
imshow(frame_unit1_offset.cdata);
T=title(unit1_text,'FontSize',14,'FontName','Arial','FontWeight','bold','Color',COL(2,:));
axes(ha(5));
imshow(frame_unit1_slope.cdata);

for i=1:2*frate   
    if(i<=D)
        axes(Abox);        
        R.FaceColor=[0 0 0 1-d_05s(i)];
    end
    writeVideo(vidout,getframe(Fall));
end

%% unit 2 movie
Abox=axes('Position',pos{3}.*[1 1 1 1.2],'Color','none','XColor','none','YColor','none');
R=rectangle('Position',[0 0 1 1],'LineStyle','none','FaceColor',[0 0 0 1]);
k=1;
for j=1:2
    for i=[1:(2*K-1)]+floor(K/2)
        axes(ha(3));
        frk=abs(K-1-abs(i-K))+1;        
        imshow(frame_unit2(frk).cdata);
        T=title(unit2_text,'FontSize',14,'FontName','Arial','FontWeight','bold','Color',COL(3,:));
        if(k<=D)
            axes(Abox);    
            R.FaceColor=[0 0 0 1-d_05s(k)];
            k=k+1;
        end        
        writeVideo(vidout,getframe(Fall));
    end
end

Abox=axes('Position',pos{6}.*[1 1 1 1],'Color','none','XColor','none','YColor','none');
R=rectangle('Position',[0 0 1 1],'LineStyle','none','FaceColor',[0 0 0 1]);
axes(ha(3));
ha(3).clo;
imshow(frame_unit2_offset.cdata);
T=title(unit2_text,'FontSize',14,'FontName','Arial','FontWeight','bold','Color',COL(3,:));
axes(ha(6));
imshow(frame_unit2_slope.cdata);

for i=1:2*frate   
    if(i<=D)
        axes(Abox);        
        R.FaceColor=[0 0 0 1-d_05s(i)];
    end
    writeVideo(vidout,getframe(Fall));
end

close(vidout);