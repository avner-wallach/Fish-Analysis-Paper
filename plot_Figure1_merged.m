% make panels for new Figure 1: Afferent input in freely moving fish
%% methodology panels' data
adir='20190624'; seg=6; fnum=9; ss=13760; spkgroup=1; chid=1; rastind=[1:2];
% script for visualizing freely moving recordings
apath='Z:\analysis_data';
if(~exist('eod2'))
    load([apath,filesep,adir,filesep,'output'],'ops','eod','file');
    eod2=eod;
    ops2=ops;
    file2=file;
    load([apath,filesep,adir,filesep,adir,'_archieve'],'frame');
end
path(path,'Z:\GitHub\Fish-Model');
path(path,'Z:\GitHub\Fish-Pipeline')
if(~exist('eod1'))
    load('Z:\analysis_data\20190131\output.mat')
    eod1=eod;
    ops1=ops;
    file1=file;    
end
load('Z:\mormyrid_data\fish_images\fish_silhouette.mat');
str_image=IMGS;
setenv('FRAMESCALE',num2str(ops2.framescale));
setenv('SESSDATE','20190131');
setenv('DATAPATH','Z:\mormyrid_data');

setenv('PXLSIZE',num2str(ops2.pxlsize)); %pxl to mm conversion
pxl2mm=ops2.pxlsize;
fdate=num2str(ops2.seg(seg).dates(1));
dname=[ops2.datapath,'\',fdate,'\'];
fname=ops2.seg(seg).files(fnum); %file to take
cvar=ops2.seg(seg).LFPgroups{ops2.seg(seg).Spkgroups(spkgroup)};
datafile=[dname,'data_',num2str(fname)];
framenum=800; %number of frames
frameind=ss+[0:framenum];

%% get video and data files
% if(~exist('data'))
%     load(datafile);        %for traces
% end
% 
% if(~exist('amp'))
%     chns=find(ops2.outchans);
%     ch=chns(cvar);
%     sesspath=[ops2.datapath,'\',fdate,'\'];
%     [t,amp]=read_bonsai_binary([sesspath,'amp_',num2str(fname)],ops2.samplerate,ops2.chan_num,ops2.blocksize,ch,'adc');
%     t=t-data.FILE.offset;
%     a=amp(:,1);
% end
%% start figure
% bgc=[0 0 0];
bgc=[1 1 1];
fontsize=7;
nline=0.25;
wline=.5;
msize=3;

Fall=figure;
set(Fall,'Units','centimeters');
set(Fall,'OuterPosition', [2  2  17.2  12],'Color',bgc);

[ha, pos] = tight_subplot(2, 1, [.1 .0],[.05 .0],[.0 .0],[.5 .5]); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]

%% first row
COL=colormap(ha(1),'lines');
% COL(4,:)=[236 157 118]/255; %afferent color
COL(4,:)=[226 139 138]/255; %afferent color
COL(5,:)=[64 128 0]/255;%[57 181 74]/255; output color
COL(6,:)=[210 148 37]/255;%COL(7,:);%interneuron
COL(7,:)=[86 124 141]/255;%output color

%% first row
[h2,p2]=tight_subplot(1, 5, [.0 .04],[.0 .0],[.00 .02],[],[.28 .25 .1 .1 .27],ha(1)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]

%% get snapshot images
getimages=0;
if(getimages)
    [r,th]=convert_exo_polar(frame(2).data(:,1),frame(2).data(:,2),...
        frame(2).data(:,3),file1(2).circle,pxl2mm);
    % ind=find(inrange(r,[49 51]) & inrange(th,pi/2+[-pi/20 pi/20]));
    x=r.*sin(th); y=r.*cos(th);
    t=frame(2).data(:,9);
    for i=1:3
%         figure(Fre);
%         h=impoly();
%         p=h.getPosition;
%         ind = find(inpolygon(x,y,p(:,1),p(:,2)));
%         %get indices for frames
%         [fileidx{i},frameidx{i}]=get_frame_indices(ind,frame(2).ind0);
%         T=mat2cell(t(ind),cellfun(@numel,frameidx{i}));
        [m,im]=max(cellfun(@numel,frameidx{i}));
        figure;
        [frames{i},ax_grid] = get_fish_image(ops1.seg(2).files(fileidx{i}(im)),frameidx{i}{im},'headfix','off');
    end
    save('fish_frames','frames','fileidx','frameidx','-v7.3');
else
    if(~exist('frames'));
        load('fish_frames','frames','fileidx','frameidx','p');
    end
end

%% get joint figure of all fish
F(:,:,:,1)=frames{1}(9).cdata;
F(:,:,:,2)=frames{2}(10).cdata;
F(:,:,:,3)=frames{3}(14).cdata;
D=double(F);
th=45;
M=median(D,4);
x=1:size(M,2); y=1:size(M,1);
[X,Y]=meshgrid(x,y);
x0=406; y0=406; r=size(M,1)/2;
B0=(((X-x0).^2+(Y-y0).^2)<r.^2);
B1=(abs(D(:,:,:,1)-M)>th);
B2=(abs(D(:,:,:,2)-M)>th);
B3=(abs(D(:,:,:,3)-M)>th);
c1(:,:,1)=57/255*ones(size(B1(:,:,1))); c1(:,:,2)=116/255*ones(size(B1(:,:,1))); c1(:,:,3)=254/255*ones(size(B1(:,:,1)));
% c2(:,:,1)=64/255*ones(size(B1(:,:,1))); c2(:,:,2)=202/255*ones(size(B1(:,:,1))); c2(:,:,3)=141/255*ones(size(B1(:,:,1)));
c2(:,:,1)=168/255*ones(size(B1(:,:,1))); c2(:,:,2)=205/255*ones(size(B1(:,:,1))); c2(:,:,3)=57/255*ones(size(B1(:,:,1)));
c3(:,:,1)=217/255*ones(size(B1(:,:,1))); c3(:,:,2)=189/255*ones(size(B1(:,:,1))); c3(:,:,3)=39/255*ones(size(B1(:,:,1)));
FF1=(255-D(:,:,:,1)).*c1;
FF2=(255-D(:,:,:,2)).*c2;
FF3=(255-D(:,:,:,3)).*c3;
figure;
%axes(h2(3));
I=image(uint8(M));
set(I,'AlphaData',B0);
hold on;
I=image(uint8(FF1));
set(I,'AlphaData',B1(:,:,1).*B0);
I=image(uint8(FF2));
set(I,'AlphaData',B2(:,:,1).*B0);
I=image(uint8(FF3));
set(I,'AlphaData',B3(:,:,1).*B0);
set(gca,'XColor','none','YColor','none');
axis('image');
axes(h2(1));
h2(1).clo;
I=image(uint8(M));
set(I,'AlphaData',B0);
hold on;
I=image(uint8(FF1));
set(I,'AlphaData',B1(:,:,1).*B0);
I=image(uint8(FF2));
set(I,'AlphaData',B2(:,:,1).*B0);
I=image(uint8(FF3));
set(I,'AlphaData',B3(:,:,1).*B0);
set(gca,'XColor','none','YColor','none');
axis('image');
% h2(1).Visible='off';
% 
%% plot ex afference, re afference maps
clim=[-1 1];
axes(h2(5));
h2(5).clo;
% lfp_cmap=ones(4,1)*COL(4,:);
[Mz,N0,s1]=plot_tuning_polar(eod1(2),file1(2),18,'t_col',6,'clim',clim,...
    'image',str_image.cdata,'fontsize',fontsize,'mfunc','slope','ops',ops1,'bgcol',bgc,...
    'upsamp',0.25,'r_nbins',10,'th_nbins',15,'poly',p); %20-20?
set(h2(5),'ALim',[1 2.5]);

[hlfp,plfp]=tight_subplot(3,2,[.01 .0],[0 0],[0 0],[],[.35 .65],h2(2));
I=[1 4 2 5 3 6];
for i=1:numel(I)
    hlfp(i).clo;
%     set(ca(i),'XLim',[-10 40]);    
    copyaxes(s1.a_out(I(i)),hlfp(i),1);
    set(hlfp(i),'Colormap',s1.a_out(I(i)).Colormap,'Color','none');
end
h2(3).clo;
copyaxes(s1.a_out(end),h2(3),1);
set(h2(3),'YTick',[-1 -.5 0 .5 1],'YTickLabelMode','auto');

% hold on;
axes(h2(4));
m=Mz(:);
n=N0(:);
n=n/max(n)*50+1;
S=scatter(randn(size(m))*.1,m,25,m,'filled','MarkerEdgeColor',0.25*[1 1 1],'MarkerFaceAlpha',1,'LineWidth',.1);
set(h2(4),'Xlim',[-1 1],'Ylim',[-1 1],'YTickLabel',[],'Clim',[-1 1],'XAxisLocation','origin','XTick',[]);
set(h2(4),'Colormap',h2(5).Colormap,'Color','none');

%% third row
% [h3,p3]=tight_subplot(1, 3, [.0 .05],[.0 .00],[.00 .0],[],[7/16 7/16 1/8],ha(3)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]
[h3,p3]=tight_subplot(1, 5, [.0 .04],[.0 .0],[.00 .02],[],[.28 .25 .1 .1 .27],ha(2)); %Nh, Nw, [gap_h gap_w], [lower upper], [left right]

h33=axes('Position',p3{1});
[hmod,pmod]=tight_subplot(3,2,[.01 .0],[0 0],[0 0],[],[.35 .65],h3(2));
set(hmod([1 3 5]),'Visible','off');
h3(1).Position(1)=p3{1}(1)+pmod{1}(3)+.05;
h3(1).Position(3)=p3{1}(3)-.05;
h33.Position(1)=p3{1}(1)+.05;
h33.Position(3)=pmod{1}(3);
%% reafference variance
load('Z:\analysis_data\ELL__Models_Database.mat');
COL2=brewermap(12,'Paired');
% reafference non uniformity
axes(h33);
% h3(3).Position=p3{3}+[0.0275 0.0125 -0.0275 -0.0375];
h33.clo;
% CC=COL2([2 1 8 7],:);
CC(1,:)=.1*[1 1 1];
CC(2,:)=.7*[1 1 1];
% CC=[COL(4,:); lighter(COL(4,:),.25)];
H=plot([1;2]*ones(size(groups.Tail_Var')),[groups.Tail_Var'; groups.Tail_Var_Ctl'],'-');
set(H(:),'MarkerSize',msize,'MarkerEdgeColor',CC(1,:),'MarkerFaceColor','none','Color',CC(2,:));
hold on;
H=scatter([ones(size(groups.Tail_Var));2*ones(size(groups.Tail_Var))],[groups.Tail_Var; groups.Tail_Var_Ctl],msize,CC(1,:),'filled');
set(H(:),'MarkerEdgeColor',CC(1,:),'MarkerFaceColor',CC(1,:),'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
set(h33,'FontSize',fontsize,'Xlim',[0.5 2.5],'Ylim',[0 0.12],'YTick',[0 .1 .2],'YTickLabel',{'0','0.1','0.2'},'XTick',[],'Color','none');
set(h33,'TickLength',[0.025 0.],'TickDir','both','Box','off');


%%
ftemp=figure; atmp=axes;
[FF,Z,T,R,TH,N,P1,P2]=get_tail_wall([atmp h3(5)],.75*pi,1);
close(ftemp);
h3(5).CLim=clim;
colorbar('Location','east','AxisLocation','out','LineWidth',nline,'TickLength',0.025,'TickDirection','both');
CMAP=colormap(h3(5));
CN=size(CMAP,1);
m=P1(:);
n=N(:);
n=n/max(n)*5+1;
row{1}=[3]; col{1}=[12]; %2 10
row{2}=[2]; col{2}=[17];
row{3}=[12]; col{3}=[9];

for i=1:3
    z{i}=[];
    for j=1:numel(row{i})
        z{i}(j,:)=permute(Z(col{i}(j),row{i}(j),:),[3 1 2]);
    end
end

T=T/(max(T(:)));
for i=1:3
    axes(hmod(i*2));
    hmod(i*2).clo;
    TT=repmat(T',size(z{i},1),1);
    tt=TT(~isnan(z{i})); zz=z{i}(~isnan(z{i}));
    zz=zz+0.1*randn(size(zz));
    S=scatter(tt,zz,1.677,tt,'filled');
    set(S,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',1,'MarkerFaceColor',0.5*[1 1 1]);    
    hold on;
    f=fit(tt(:),zz(:),'poly1');
    sl1(i)=f.p1;
    cval=min(max((f.p1-clim(1))/diff(clim),0),1);
    MCOL(i,:)=interp1(linspace(0,1,CN),CMAP,cval);
    pp = predint(f,T,0.95,'functional','on');
    A1=patch([T' fliplr(T')],[pp(:,2)' fliplr(pp(:,1)')],MCOL(i,:),'LineStyle','none','FaceColor',MCOL(i,:),'FaceAlpha',0.5);
    hold on;
    H=plot(T,f(T)); set(H,'Color',MCOL(i,:),'LineWidth',3*wline); legend('off');
    H=plot(T,pp(:,1)); set(H,'Color',0.5*[1 1 1],'LineWidth',nline); legend('off');
    H=plot(T,pp(:,2)); set(H,'Color',0.5*[1 1 1],'LineWidth',nline); legend('off');
    
%     XX=linspace(min(T),max(T),1e2);
%     hold on;
%     H=plot(XX,f(XX)); set(H,'Color',MCOL(i,:),'LineWidth',wline); legend('off');
    
    legend('off');    
    xlabel('');
    ylabel('');
    
    set(gca,'XAxisLocation','origin','YAxisLocation','origin',...
        'Ylim',[-1 1],'Xlim',[-2 2],'YTick',[],'XTick',[]);    
%     set(gca,'Color',bgcol,'YColor',1-bgcol,'XColor',1-bgcol,'Colormap',CC);    
end
axes(h3(3));
h3(3).clo;
for i=1:3
    H=plot(i,sl1(i),'o');
    set(H,'MarkerSize',5,'MarkerFaceColor',MCOL(i,:),'MarkerEdgeColor',0*[1 1 1],'LineWidth',1);
    hold on;
end
  set(h3(3),'Color',bgc,'YColor',1-bgc,'XColor',1-bgc,'Xlim',[.5 3.5],'Ylim',[-1 1],'XTick',[],'YTick',[-1 -.5 0 .5 1],'XAxisLocation','origin','Box','off');    

axes(h3(4));
h3(4).clo;
S=scatter(randn(size(m))*.1,m,25,m,'filled','MarkerEdgeColor',0*[1 1 1],'MarkerFaceAlpha',1,'LineWidth',.1);
set(h3(4),'Xlim',[-1 1],'Ylim',[-1 1],'Clim',[-1 1],'XAxisLocation','origin','XTick',[],'YTickLabel',[],'Position',p3{4});
colormap(gca,CMAP);

C=findobj('Type','Colorbar');
set(C(:),'Location','east','AxisLocation','out','LineWidth',nline,'TickLength',0.025,'TickDirection','both','FontSize',fontsize,'FontWeight','normal');
for i=1:numel(C)
    C(i).Position(3)=C(i).Position(3)/2;
%     C(i).Position(4)=C(i).Position(4)*.7;
end

axes(h3(1));
h3(1).clo;
plot_map('tail_angle',.225*pi,'tank_radius',23,'wall_dist',5,'wall_angle',.95*pi,...
    'p_density',20,'fish_length',15,'grid_M',500,'r_max',50,'tail_p',0.475,...
    'plot_potential',1,'plot_field',1,'plot_lfield',0,'reflection',1,'mpos','.','mneg','.','axes',h3(1));
% axis('image');
set(gca,'Xlim',[-20 15],'YLim',[-28 8]);
%% rearange
h2(1).Color='none';
h3(4).Color='none';
h3(1).Color='none';
h2(3).Position(1)=p2{3}(1)+.035;
h2(3).Color='none';
h2(5).Position(1)=p2{5}(1)-.035;
h2(5).Color='none';

h3(3).Position(1)=p3{3}(1)+.035;
h3(3).Color='none';
h3(5).Position(1)=p3{5}(1)-.035;
h3(5).Color='none';

%% remove images before copy
% hs={h2(2),h31(1),h31(2),h33(1),h33(2)};
srf=findobj('Type','Surface');
ptch=findobj('Type','Patch');
img=findobj('Type','Image');
set(srf,'Visible','off');
set(img,'Visible','off');
for i=1:numel(ptch)
    if(ptch(i).Parent==h3(5))
       set(ptch(i),'Visible','off');
    end
end
set(Fall,'Color','none');
print('-clipboard','-dmeta');
%% remove all but figures before copy to bmp
hall=findobj();
for i=1:numel(hall)
    if(~strcmp(hall(i).Type,'root') & ~strcmp(hall(i).Type,'figure'))
        set(hall(i),'Visible','off');
    end
end
srf=findobj('Type','Surface');
ptch=findobj('Type','Patch');
img=findobj('Type','Image');
set(srf,'Visible','on');
set(img,'Visible','on');
for i=1:numel(ptch)
    if(ptch(i).Parent==h3(5))
       set(ptch(i),'Visible','on');
    end
end
print('-clipboard','-dbitmap','-r720'); %copy maps

        
