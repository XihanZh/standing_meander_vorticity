%================fig1======================
load manuscript.mat
load cooreal.mat
%-thermal colormap is downloaded from 
%https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
%--upload m files in github_repo
% load gfdl.mat
% d=rdmds('Depth');
close all
x0=10;
y0=10;
width=800;
height=3200;
set(gcf,'position',[x0,y0,width,height])
% fig=gcf;
ax(1)=subplot(3,1,1);
imagesc(x(46:end-45),y(26:end-25),d(46:end-45,26:end-25)'),axis xy
colormap(ax(1),parula)
colorbar
caxis([0 6000])
colorbar('Ticks',0:1000:6000)
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
% ylabel('Latitude','fontsize',14)
tobj=title('(a)','fontsize',14);
% tobj.Position=[x(1)+0.4,y(end)+0.1];
tobj.Position=[x(46)+0.4,y(end-25)+0.1];
set(gca,'fontsize',14)
ylabel(colorbar,'Bathymetry [m]')

ax(2)=subplot(3,1,2);
Ts(Ts==0)=nan;
imagesc(x(46:end-45),y(26:end-1),Ts(46:end-45,26:end-1)'),axis xy
colorbar
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
cmap=cmocean('thermal');
colormap(ax(2),cmap)
hold on
a=15;
xnew=x(46:a:end-45);
ynew=y(26:a:end-25);
[X,Y]=meshgrid(xnew,ynew);
uq=us(46:a:end-45,26:a:end-25)';
vq=vs(46:a:end-45,26:a:end-25)';
[X1,Y1]=meshgrid(xnew,y(end-15));
%Then a reference arrow velocity is defined at x=0.2 and y=2.2, while rest of points are left as zero:
u_rf=zeros(1,length(xnew));
v_rf=zeros(1,length(xnew));
u_rf(1,2)=2;

quiver([X;X1],[Y;Y1],[uq;u_rf],[vq;v_rf],'color',[0 0 0],'LineWidth',1)
% ylabel('Latitude','fontsize',14)
set(gca,'fontsize',14)
text(x(47),y(end-25),'Arrow scale: 2m/s')
ylabel(colorbar,'SST [^{o}C]')

set(gca,'fontsize',14)
tobj=title('(b)','fontsize',14);
tobj.Position=[x(46)+0.4,y(end-1)+0.1];


ax(3)=subplot(3,1,3);
imagesc(x(46:end-45),y(26:end-25),eke(46:end-45,26:end-25)'),axis xy
colorbar %caxis([-5e-4 5e-4])%caxis([-8e-3 8e-3])
set(gca,'fontsize',14)
caxis([0 0.20])
colormap(ax(3),'hot')
set(ax(3),'xtick',[],'ytick',[]);
ylabel(colorbar,'EKE [m^{2}s^{-2}]')
%colormap(h1,map)
v=-0.5:0.1:0.5;

[X,Y]=meshgrid(x(46:end-45),y(26:end-25));
h2=axes('position',get(ax(3),'position'),'color','none','fontsize',14); %,'XLim',[145 175],'YLim',[48 60],'YDir', 'reverse');
hold on
[cn,hn]=contour(X,Y,etabar(46:end-45,26:end-25)',v,'showtext','on','LineWidth',2,'parent',h2);
clabel(cn,hn,-0.5:0.2:0.5,'color','c','fontsize',10)
colormap(h2,[1,1,1])
B2=[230,460,125,343];%large meander
rr=B2;
plot(x(rr(1))*ones(1,length(rr(3):rr(4))),y(rr(3):rr(4)),'b','LineWidth',2)
plot(x(rr(2))*ones(1,length(rr(3):rr(4))),y(rr(3):rr(4)),'b','LineWidth',2)
plot(x(rr(1):rr(2)),y(rr(3))*ones(1,length(rr(1):rr(2))),'b','LineWidth',2)
plot(x(rr(1):rr(2)),y(rr(4))*ones(1,length(rr(1):rr(2))),'b','LineWidth',2)
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
tobj=title('(c)','fontsize',14);
tobj.Position=[x(46)+0.4,y(end-25)+0.1];
% xlabel('Longitude ','fontsize',14)
% ylabel('Latitude','fontsize',14)
set(gca,'fontsize',14)
print -dpng fig1_depthandEKE.png
%===========fig2=======================
load manuscriptvort.mat
load cooreal.mat
load blue_red_saturated.mat
load manuscript.mat
x0=10;
y0=10;
width=800;
height=3200;
set(gcf,'position',[x0,y0,width,height])

ax(1)=subplot(3,1,1);
imagesc(x(45:end-44),y(26:end-25),1e11*sadv(44:end-44,25:end-25)'),axis xy
colorbar('Ticks',-4:2:4)
colormap(ax(1),map)
% ylabel(colorbar,'\boldmath$\bf{u}\cdot\nabla\zeta$','Interpreter','latex','fontsize',14)%'Rotation',0)
ylabel(colorbar,'$\bf{u}\cdot\nabla\zeta$','Interpreter','latex','fontsize',18)%'Rotation',0)
% caxis([-5e-11 5e-11])
caxis([-5 5])
set(ax(1),'XTick',[],'YTick',[])
h1=axes('position',get(ax(1),'position'),'color','none','fontsize',14);
hold on
contour(x(45:end-44),y(26:end-25),etabar(45:end-44,26:end-25)','ShowText','on')
colormap(h1,[0,0,0])
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
tobj=title('(a)','fontsize',14);
tobj.Position=[x(46)+0.4,y(end-25)+0.1];
% ylabel('Latitude','fontsize',14)
set(gca,'fontsize',14)

ax(2)=subplot(3,1,2);
imagesc(x(45:end-44),y(26:end-25),1e11*stretch(44:end-44,25:end-25)'),axis xy
colorbar('Ticks',-4:2:4)
caxis([-5 5])
colormap(ax(2),map)
% ylabel(colorbar,'\boldmath$f\frac{\partial w}{\partial z}$','Interpreter','latex','fontsize',14)%'Rotation',0)
ylabel(colorbar,'$-f{\frac{\partial w}{\partial z}}$','Interpreter','latex','fontsize',20)%'Rotation',0)
set(ax(2),'XTick',[],'YTick',[])
h2=axes('position',get(ax(2),'position'),'color','none','fontsize',14);
hold on
contour(x(45:end-44),y(26:end-25),etabar(45:end-44,26:end-25)','ShowText','on')
colormap(h2,[0,0,0])
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
tobj=title('(b)','fontsize',14);
tobj.Position=[x(46)+0.4,y(end-25)+0.1];
% ylabel('Latitude','fontsize',14)
set(gca,'fontsize',14)


ax(3)=subplot(3,1,3);
imagesc(x(45:end-44),y(26:end-25),1e11*bv(44:end-44,25:end-25)'),axis xy
colorbar('Ticks',-0.4:0.2:0.4)
caxis([-0.5 0.5])
colormap(ax(3),map)
% ylabel(colorbar,'\boldmath$\beta v$ $(\times10^{11}s^{-2})$','Interpreter','latex','fontsize',14,'fontweight','bold')%'Rotation',0)
ylabel(colorbar,'$\beta v$','Interpreter','latex','fontsize',20)%'Rotation',0)
set(ax(3),'XTick',[],'YTick',[])
h3=axes('position',get(ax(3),'position'),'color','none','fontsize',14);
hold on
contour(x(45:end-44),y(26:end-25),etabar(45:end-44,26:end-25)','ShowText','on')
colormap(h3,[0,0,0])
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
tobj=title('(c)','fontsize',14);
tobj.Position=[x(46)+0.4,y(end-25)+0.1];
% ylabel('Latitude','fontsize',14)
% xlabel('Longitude','fontsize',14)
set(gca,'fontsize',14)
% colormap(map)
print -dpng fig2_maps.png
%=====figure3==============
load manuscriptvort.mat
load cooreal.mat
load blue_red_saturated.mat
load vortsec54S.mat
% d=rdmds('Depth');
% la=233;
la=200;
[X,Z]=meshgrid(x(45:end-45),zc(1:168));
close all
x0=10;
y0=10;
width=800;
height=3200;
set(gcf,'position',[x0,y0,width,height])

subplot(3,1,1)
pcolor(X,Z,1e11*sadvsec(45:end-44,1:168)')
shading flat
caxis([-5 5])
colorbar('Ticks',-4:2:4)
% ylabel(colorbar,'\boldmath$\bf{u}\cdot\nabla\zeta$','Interpreter','latex','fontsize',14)%'Rotation',0)
ylabel(colorbar,'$\bf{u}\cdot\nabla\zeta$','Interpreter','latex','fontsize',18)%'Rotation',0)
hold on
plot(x(45:end-45),-d(45:end-45,la),'k','LineWidth',1)
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
ylabel('Depth [m]','fontsize',14)
set(gca,'fontsize',14)
set(gca,'Layer','top')
tobj=title('(a)','fontsize',12);
tobj.Position=[x(45)+0.4,zc(1)+0.1];

subplot(3,1,2)
pcolor(X,Z,1e11*stretchsec(45:end-44,1:168)')
shading flat
caxis([-5 5])
colormap(map)
colorbar('Ticks',-4:2:4)
% ylabel(colorbar,'\boldmath$f\frac{\partial w}{\partial z}$','Interpreter','latex','fontsize',14)%'Rotation',0)
ylabel(colorbar,'$-f{\frac{\partial w}{\partial z}}$','Interpreter','latex','fontsize',20)%'Rotation',0)
hold on
plot(x(45:end-45),-d(45:end-45,la),'k','LineWidth',1)
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
ylabel('Depth [m]','fontsize',14)
set(gca,'fontsize',14)
set(gca,'Layer','top')
tobj=title('(b)','fontsize',12);
tobj.Position=[x(45)+0.4,zc(1)+0.1];

subplot(3,1,3)
pcolor(X,Z,1e11*bvsec(45:end-44,1:168)')
shading flat
caxis([-0.5 0.5])
colormap(map)
colorbar('Ticks',-0.4:0.2:0.4)
% ylabel(colorbar,'\boldmath$\beta v$','Interpreter','latex','fontsize',14,'fontweight','bold')%'Rotation',0)
ylabel(colorbar,'$\beta v$','Interpreter','latex','fontsize',20)%'Rotation',0)
hold on
plot(x(45:end-45),-d(45:end-45,la),'k','LineWidth',1)
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
ylabel('Depth [m]','fontsize',14)
set(gca,'fontsize',14)
set(gca,'Layer','top')
tobj=title('(c)','fontsize',12);
tobj.Position=[x(45)+0.4,zc(1)+0.1];
colormap(map)
print -dpng fig3_verticalsections.png
%-----figure4
load manuscriptvort_ideal.mat
load blue_red_saturated.mat
x=x*1e-3;
y=y*1e-3;
la=116;
[X,Y]=meshgrid(x(2:end-1),y(2:end-1));
close all
x0=10;
y0=10;
width=800;
height=3200;
set(gcf,'position',[x0,y0,width,height])
ax(1)=subplot(4,1,1);
imagesc(x(2:end-1),y(2:end-1),d(2:end-1,2:end-1)'),axis xy
colorbar
caxis([2000 3000])
colorbar('Ticks',2000:200:3000)
ylabel(colorbar,'Bathymetry [m]','fontsize',12)
% colormap(ax(1),'hot')
set(ax(1),'xtick',[],'ytick',[]);
h2=axes('position',get(ax(1),'position'),'color','none','fontsize',14); 
hold on
contour(X,Y,etabar(2:end-1,2:end-1)','LineWidth',2,'ShowText','on','parent',h2)
colormap(h2,[0,0,0])
ylabel('Distance [km]','fontsize',14)
% rm=[95,170,35,170];
rm=[77,163,38,190];
plot(x(rm(1):rm(2)),y(rm(3))*ones(1,length(rm(1):rm(2))),'r','LineWidth',2)
plot(x(rm(1):rm(2)),y(rm(4))*ones(1,length(rm(1):rm(2))),'r','LineWidth',2)
plot(x(rm(1))*ones(1,length(rm(3):rm(4))),y(rm(3):rm(4)),'r','LineWidth',2)
plot(x(rm(2))*ones(1,length(rm(3):rm(4))),y(rm(3):rm(4)),'r','LineWidth',2)
xticks(0:200:2000)
set(gca,'fontsize',14)
title('(a)','position',[x(2)+30 y(end-1)+20])

ax(2)=subplot(4,1,2);
% x1=[x,2*x(end)-x(end-1)];
% zc1=[-0.1;zc];
[X,Z]=meshgrid(x,zc);
advnew=nan(size(adv));
advnew(1:end-1,1:end-1)=adv(2:end,2:end);
pcolor(X,Z,1e11*advnew')
% pcolor(X,Z,1e11*adv')
shading flat
% colorbar
% colorbar('Ticks',-3:1:3)
colormap(ax(2),map)
% caxis([-1e-11 1e-11])
caxis([-1 1])
hold on
plot(x,-d(:,la),'k','LineWidth',2)
set(gca,'Layer','top','fontsize',14)
xticks(0:200:2000)
% xticklabels({'135^{\circ}E','140^{\circ}E','145^{\circ}E','150^{\circ}E',...
%     '155^{\circ}E','160^{\circ}E','165^{\circ}E'})
title('(b)','position',[x(2)+30 zc(1)])
ylabel('Depth [m]')
% ylabel(colorbar,'\boldmath$\bf{u}\cdot\nabla\zeta$','Interpreter','latex','fontsize',14)%'Rotation',0)
ylabel(colorbar,'$\bf{u}\cdot\nabla\zeta$','Interpreter','latex','fontsize',18)%'Rotation',0)
ax(3)=subplot(4,1,3);
[X,Z]=meshgrid(x,zc);
strecthnew=nan(size(stretch));
strecthnew(1:end-1,1:end-1)=stretch(2:end,2:end);
pcolor(X,Z,1e11*strecthnew')
shading flat
colorbar
colormap(ax(3),map)
caxis([-1 1])
hold on
plot(x,-d(:,la),'k','LineWidth',2)
set(gca,'Layer','top','fontsize',14)
title('(c)','position',[x(2)+30 zc(1)])
ylabel('Depth [m]')
% ylabel(colorbar,'\boldmath$f\frac{\partial w}{\partial z}$','Interpreter','latex','fontsize',14)%'Rotation',0)
ylabel(colorbar,'$-f{\frac{\partial w}{\partial z}}$','Interpreter','latex','fontsize',20)%'Rotation',0)


ax(4)=subplot(4,1,4);
[X,Z]=meshgrid(x,zc);
bvnew=nan(size(bv));
bvnew(1:end-1,1:end-1)=bv(2:end,2:end);
pcolor(X,Z,1e11*bvnew')
shading flat
colorbar
colormap(ax(4),map)
caxis([-1 1])
hold on
plot(x,-d(:,la),'k','LineWidth',2)
set(gca,'Layer','top','fontsize',14)
title('(d)','position',[x(2)+30 zc(1)])
ylabel('Depth [m]')
xlabel('Distance [km]')
% ylabel(colorbar,'\boldmath$\beta v$','Interpreter','latex','fontsize',14,'fontweight','bold')%'Rotation',0)
ylabel(colorbar,'$\beta v$','Interpreter','latex','fontsize',20)%'Rotation',0)
print -dpng fig4_idealised.png

%==figure5====================
H=4e3;
dz=20;
z=-H:dz:0;
nz=length(z);
u0=0.01;
U1=0.2;
U2=0.01;
ubar=zeros(1,nz);
ubar(1:150)=U2;
ubar(151:end)=U1;
ub=0.01;
us=0.3;
sigma=1e3;
% u=(ut*exp(z/sigma))/(exp(-H/sigma))+ub;
u=(us-ub)*(1-exp((z+H)/sigma))/(1-exp(H/sigma))+ub;
close all
figure
plot(u0*ones(1,nz),z,'k','LineWidth',2)
hold on
plot(ubar,z,'r--','LineWidth',2)
plot(u,z,'LineWidth',2)
grid on
ylabel('Depth [m]','fontsize',14)
xlabel('U [m/s]','fontsize',14)
legend({'uniform','step-function','continuous'},'fontsize',14,'location','best')
set(gca,'fontsize',14)
print -dpng fig5_QGprofiles.png

%=====figure6=============
load blue_red_saturated.mat
N=1e-3;
f=-1e-4;
H=4e3;
beta=1e-11;
ubar=0.01;
dz=20;
z=-H:dz:0;
h0=200;
L=3e4;
k=2*pi/L;
lamsq=(N^2/f^2)*(beta/ubar-k^2);
caplam=abs(sqrt(lamsq));
phi=((N^2)*h0)/(f*caplam)*cosh(caplam*z)/sinh(caplam*H); % analytical solution
nx=201;
% L1=5e5;
L1=3e5;
x=linspace(0,L1,nx);
[X,Z]=meshgrid(x,z);
psi=(((N^2)*h0)/(f*caplam)*cosh(caplam*Z)/sinh(caplam*H)).*exp(1i*k*X);
U=ubar*ones(1,length(z));
v=-k*repmat(phi',[1,length(x)]).*sin(k*X);
adv=-k^2*repmat(U',[1,length(x)]).*v;
betav_s=-k*beta*phi;
adv_s=k^3*ubar*phi;
phiz=(phi(1:end-1)-phi(2:end))/dz;
phizz=(phiz(1:end-1)-phiz(2:end))/dz;
stretch_s=-(ubar*k*f^2/N^2)*phizz;
close all
subplot(2,1,1)
pcolor(X*1e-3,Z,adv)
shading flat
colorbar
% caxis([-5e-10 5e-10])
caxis([-1e-10 1e-10])
colormap(map1)
% ylabel(colorbar,'\boldmath$\bf{u}\cdot\nabla\zeta$','Interpreter','latex','fontsize',14)
ylabel(colorbar,'$\bf{u}\cdot\nabla\zeta$','Interpreter','latex','fontsize',18)%'Rotation',0)
tobj=title('(a)','fontsize',12);
tobj.Position=[x(1)+10,z(end)+50];
ylabel('Depth [m]','fontsize',14)
xlabel('Distance [km]','fontsize',14)
set(gca,'Layer','top','fontsize',14)

subplot(2,1,2)
plot(betav_s,z,'LineWidth',2)
hold on
plot(adv_s,z,'LineWidth',2)
plot(stretch_s,z(2:end-1),'LineWidth',2)
grid on
% plot(betav_s(2:end-1)+adv_s(2:end-1)+stretch_s,z(2:end-1),'k','LineWidth',2)
ylabel('Depth [m]','fontsize',14)
tobj=title('(b)','fontsize',12);
tobj.Position=[-9.5e-11,0.12];
legend({'\beta v','advection','stretching'},'fontsize',14,'Location','best')
xlabel('Vorticity terms [s^{-2}]','fontsize',14)
set(gca,'fontsize',14)
print -dpng fig6_QGsmall.png
%-----figure 7
load blue_red_saturated.mat
N=1e-3;
f=-1e-4;
H=4e3;
beta=1e-11;
ubar=0.01;
dz=20;
z=-H:dz:0;
h0=1e3;
L=5e5;
k=2*pi/L;
lamsq=(N^2/f^2)*(beta/ubar-k^2);
caplam=abs(sqrt(lamsq));
phi=-((N^2)*h0)/(f*caplam)*cos(caplam*z)/sin(caplam*H);
nx=201;
L1=1e6;
x=linspace(0,L1,nx);
[X,Z]=meshgrid(x,z);
U=ubar*ones(1,length(z));
v=-k*repmat(phi',[1,length(x)]).*sin(k*X);
betav=beta*v;
betav_s=-k*beta*phi;
adv_s=k^3*ubar*phi;
phiz=(phi(1:end-1)-phi(2:end))/dz;
phizz=(phiz(1:end-1)-phiz(2:end))/dz;
stretch_s=-(ubar*k*f^2/N^2)*phizz;

close all
subplot(2,1,1)
% pcolor(X,Z,adv)
pcolor(1e-3*X,Z,betav)
shading flat
colorbar
% caxis([-2e-11 2e-11])
caxis([-4e-12 4e-12])
colorbar('Ticks',-4e-12:2e-12:4e-12)
colormap(map1)
ylabel(colorbar,'\boldmath$\beta v$','Interpreter','latex','fontsize',14,'fontweight','bold')%'Rotation',0)
tobj=title('(a)','fontsize',12);
tobj.Position=[x(1)+10,z(end)+50];
ylabel('Depth [m]','fontsize',14)
xlabel('Distance [km]','fontsize',14)
set(gca,'Layer','top','fontsize',14)

subplot(2,1,2)
plot(betav_s,z,'LineWidth',2)
hold on
plot(adv_s,z,'LineWidth',2)
plot(stretch_s,z(2:end-1),'LineWidth',2)
grid on
ylabel('Depth [m]','fontsize',14)
tobj=title('(b)','fontsize',12);
tobj.Position=[-4.8e-12,0.12];
xlim([-5e-12 5e-12])
xticks(-5e-12:1e-12:5e-12)
legend({'\beta v','advection','stretching'},'fontsize',14,'Location','best')
xlabel('Vorticity terms [s^{-2}]','fontsize',14)
set(gca,'fontsize',14)
print -dpng fig7_QGlarge.png

%====fig8
load blue_red_saturated.mat
N=1e-3;
h0=1e3;
f=-1e-4;
H=4e3;
beta=1e-11;
dz=20;
z=-H:dz:0;
nz=length(z);
ub=0.01;
us=0.3;
sigma=1e3;
% u=(ut*exp(z/sigma))/(exp(-H/sigma))+ub;
u=(us-ub)*(1-exp((z+H)/sigma))/(1-exp(H/sigma))+ub;
L=5e5;
k=2*pi/L;
lamsq=(N^2/f^2).*(beta./u-k^2);
D2zt = -2*diag(ones(1,nz)) + diag(ones(1,nz-1), 1) + diag(ones(1,nz-1), -1);
D2zt(1,:) = 0*D2zt(1,:);
D2zt(nz,:) = 0*D2zt(nz,:);
D2zt = D2zt/dz^2;

M=N^2/f^2*(diag(beta./u)-k^2*eye(nz));
M(1,:)=0*M(1,:);
M(nz,:)=0*M(nz,:);


E=M+D2zt;
E(1,1)=-1/dz;
E(1,2)=1/dz;
E(end,end-1)=-1/dz;
E(end,end)=1/dz;

R=zeros(1,nz);
R(end)=0;
R(1)=-(N^2/f)*h0;

invE=eye(nz)/E;
phi=invE*R';
nx=201;
L1=1e6;
x=linspace(0,L1,nx);
[X,Z]=meshgrid(x,z);
v1=-k*repmat(phi,[1,length(x)]).*sin(k*X);
advreal1=-k^2*repmat(u',[1,length(x)]).*v1;


% [X1,Z1]=meshgrid(x,z(2:end-1));
v=-k*phi;
betavreal=beta*v;
advreal=-k^2*(u').*v;

phiz=(phi(1:end-1)-phi(2:end))/dz;
phizz=(phiz(1:end-1)-phiz(2:end))/dz;
% psixzz=-repmat(phizz,[1,length(x)]).*sin(k*X1);%99*201

% stretchingreal=k*(f^2/N^2)*repmat(u(2:end-1)',[1,length(x)]).*psixzz;
stretchingreal=-k*(f^2/N^2)*(u(2:end-1)').*phizz;

close all
subplot(2,1,1)
pcolor(1e-3*X,Z,advreal1)
shading flat
colorbar
% caxis([-2e-11 2e-11])
caxis([-8e-11 8e-11])
colorbar('Ticks',-8e-11:2e-11:8e-11)
colormap(map1)
% ylabel(colorbar,'\boldmath$\bf{u}\cdot\nabla\zeta$','Interpreter','latex','fontsize',14)
ylabel(colorbar,'$\bf{u}\cdot\nabla\zeta$','Interpreter','latex','fontsize',18)%'Rotation',0)
tobj=title('(a)','fontsize',12);
tobj.Position=[x(1)+10,z(end)+50];
ylabel('Depth [m]','fontsize',14)
xlabel('Distance [km]','fontsize',14)
set(gca,'Layer','top','fontsize',14)

subplot(2,1,2)
plot(betavreal,z,'LineWidth',2)
hold on
plot(advreal,z,'LineWidth',2)
plot(stretchingreal,z(2:end-1),'LineWidth',2)
grid on
ylabel('Depth [m]','fontsize',14)
tobj=title('(b)','fontsize',12);
tobj.Position=[-9.8e-11,0.12];
xlim([-1e-10 1e-10])
xticks(-1e-10:2e-11:1e-10)
legend({'\beta v','advection','stretching'},'fontsize',14,'Location','best')
set(gca,'Layer','top','fontsize',14)
xlabel('Vorticity terms [s^{-2}]','fontsize',14)
print -dpng fig8_ideallysheared.png
%-----fig9
N=1e-3;
f=-1e-4;
H=4e3;
beta=1e-11;
h0=1e3;
dz=20;
z=-H:dz:0;
nz=length(z);
U1=0.2;
U2=0.01;
ubar=zeros(1,nz);
ubar(1:150)=U2;
ubar(151:end)=U1;
L=5e5;
k=2*pi/L;
lamsq1=(N^2/f^2).*(beta/U1-k^2);
lamsq2=(N^2/f^2).*(beta/U2-k^2);
m1=sqrt(-lamsq1);
m2=abs(sqrt(-lamsq2));

H2=0.75*H;
H1=0.25*H;
alpha=(cosh(m1*H1)-U1*m1*cot(m2*H2)*sinh(m1*H1)/(U2*m2))^-1;
A=-alpha*N^2*h0/(2*f*m2*sin(m2*H2));
phianly1=2*A*cosh(m1*z(151:end));%-2000~0
C=-(N^2)*h0/(f*m2);
D=(2*A*cosh(m1*H1)+(N^2)*h0*sin(m2*H2)/(f*m2))/cos(m2*H2);
phianly2=C*sin(m2*(z(1:150)+H))+D*cos(m2*(z(1:150)+H));%-4000~-2000
phianly=[phianly2';phianly1'];



betavz=-k*beta*phianly;
advz=k^3*(ubar').*phianly;
lamdaphizz1=m1^2*phianly1;
lamdaphizz2=-m2^2*phianly2;
stretch=-(ubar').*k*(f^2/N^2).*[lamdaphizz2';lamdaphizz1'];


plot(betavz,z,'LineWidth',2)
hold on
plot(advz,z,'LineWidth',2)
plot(stretch,z,'LineWidth',2)
grid on
ylabel('Depth [m]','fontsize',14)
% tobj=title('(b)','fontsize',12);
% tobj.Position=[-1.8e-11,0.12];
xlim([-5e-11 5e-11])
xticks(-5e-11:1e-11:5e-11)
legend({'\beta v','advection','stretching'},'fontsize',14,'Location','best')
set(gca,'fontsize',14)
xlabel('Vorticity terms [s^{-2}]','fontsize',14)
print -dpng fig9_stepvort.png

%====fig10
load windsensideal.mat
close all
figure
tau=[0.5,1,2,3];
total=[mean(tothf),mean(totref),mean(totdb),mean(tottp)];
bc=[mean(bchf),mean(bcref),mean(bcdb),mean(bctp)];
bt=[mean(bthf),mean(btref),mean(btdb),mean(bttp)];
plot(tau,total,'-s','LineWidth',2)
hold on
plot(tau,bt,'-s','LineWidth',2)
plot(tau,bc,'-s','LineWidth',2)
grid on
legend({'total','barotropic','baroclinic'},'fontsize',14,'Location','best')
% xlabel('$\times\tau$','interpreter','latex','fontsize',16)
xlabel('wind multiplier','fontsize',14)
ylabel('transport [Sv]','fontsize',14)
set(gca,'fontsize',14)
print -dpng ~/Desktop/fig10_idealtrans.png

%====fig11
load windsensideal.mat
close all
figure
meander=[77,163,38,190];
rr=meander;
yyaxis left
plot([0.5,1,2,3],[mean(abs(advhf(rr(1):rr(2),rr(3):rr(4))),'all'),...
    mean(abs(advref(rr(1):rr(2),rr(3):rr(4))),'all'),mean(abs(advdb(rr(1):rr(2),rr(3):rr(4))),'all'),...
    mean(abs(advtp(rr(1):rr(2),rr(3):rr(4))),'all')],'-s','LineWidth',2)
ylabel('\boldmath$|\bf{u}\cdot\nabla\zeta|$','Interpreter','latex','fontsize',14)
hold on
yyaxis right
plot([0.5,1,2,3],[mean(lenhf,'omitnan'),mean(lenref,'omitnan'),mean(lendb,'omitnan'),mean(lentp,'omitnan')],'-s','LineWidth',2)
% ylim([2.8e6 4.8e6])
ylabel('contour length [m]','fontsize',14)
grid on
%xlabel('$\times\tau$','interpreter','latex','fontsize',14)
xlabel('wind multiplier','fontsize',14)
set(gca,'fontsize',14)
print -dpng fig11_idealmeander.png

%====fig12
load windsensideal.mat
meander=[77,163,38,190];
rr=meander;
bt=[mean(bthf),mean(btref),mean(btdb),mean(bttp)];
adv=[mean(abs(advhf(rr(1):rr(2),rr(3):rr(4))),'all'),mean(abs(advref(rr(1):rr(2),rr(3):rr(4))),'all'),...
    mean(abs(advdb(rr(1):rr(2),rr(3):rr(4))),'all'),mean(abs(advtp(rr(1):rr(2),rr(3):rr(4))),'all')];
close all
% plot(bt,1e8*adv,'s-','LineWidth',2)
plot(bt,1e8*adv,'.','MarkerSize',20)
hold on
plot(0:16,0:16,'k-')
grid on
xlabel('barotropic transport [Sv]','fontsize',14)
ylabel('\boldmath$|\bf{u}\cdot\nabla\zeta|$','Interpreter','latex','fontsize',14)
% ylabel('advection [s^{-2}]','fontsize',14)
set(gca,'fontsize',14)
xlim([0 16])
ylim([0 16])
print -dpng ~/Desktop/fig12_idealadvmeander.png

%=====fig13
load westransreal.mat
% weak=[1:5,12];
% strong=6:11;
t=1:12;
close all
yyaxis left
plot(t,trans_overtime,'s-','LineWidth',2)
hold on
% plot(t(weak),trans_overtime(weak),'k+','MarkerSize',12,'LineWidth',2)
% plot(t(strong),trans_overtime(strong),'m+','MarkerSize',12,'LineWidth',2)
plot(t,trans_bc,'s--','LineWidth',2)
ylim([110 150])
yticks(110:5:150)
% yticklabels({num2str(0),num2str(5),num2str(10),num2str(15),...
%     num2str(20),num2str(25),num2str(30),num2str(35),num2str(40)})

yyaxis right
plot(t,trans_bt,'s-','LineWidth',2)
ylim([0 40])
yticks(0:5:40)
grid on
legend({'total','baroclinic','barotropic'},'fontsize',14,'Location','bestoutside')
% legend({'total','weak','strong','BC','BT'},'fontsize',14,'Location','bestoutside')
xlabel('time [month]','fontsize',14)
title('transport [Sv]','fontsize',14)
xticks(1:12)
xticklabels({'JAN','FEB','MAR','APR',...
    'MAY','JUN','JUL','AUG','SEPT','OCT','NOV','DEC'})
xlim([1 12])
% yticks(0:30:150)
set(gca,'fontsize',14)
print -dpng fig13_current.png

%--fig14
load etawksts1month.mat
load cooreal.mat
W=[46,600-45,26,400-25];
rr=W;
v=-0.5:0.1:0.5;
close all
h1=gca;
contour(x(rr(1):rr(2)),y(rr(3):rr(4)),etabarwk(rr(1):rr(2),rr(3):rr(4))',v,'showtext','off','LineWidth',2)
colormap(h1,[0,0,0])
set(h1,'XTick',[],'YTick',[])
h2=axes('position',get(h1,'position'),'color','none','fontsize',14);
hold on
contour(x(rr(1):rr(2)),y(rr(3):rr(4)),etabarstr(rr(1):rr(2),rr(3):rr(4))',v,'showtext','off','LineWidth',2)
colormap(h2,[1,0,0])
yticks(-58:2:-50)
yticklabels({'58^{\circ}S','56^{\circ}S','54^{\circ}S',...
    '52^{\circ}S','50^{\circ}S'})
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
xlabel('Longitude','fontsize',14)
ylabel('Latitude','fontsize',14)
set(h2,'fontsize',14)
print -dpng fig14_sshsense.png

%--fig15
load vortsecsense.mat
load cooreal.mat
load blue_red_saturated.mat
d=rdmds('Depth');
la=200;
[X,Z]=meshgrid(x(45:end-45),zc(1:168));
close all
x0=10;
y0=10;
width=800;
height=4300;
set(gcf,'position',[x0,y0,width,height])

subplot(4,1,1)
pcolor(X,Z,1e11*sadvsecwkmonth0(45:end-44,1:168)')
shading flat
caxis([-5 5])
colorbar('Ticks',-4:2:4)
% ylabel(colorbar,'\boldmath$\bf{u}\cdot\nabla\zeta$','Interpreter','latex','fontsize',14)%'Rotation',0)
ylabel(colorbar,'$\bf{u}\cdot\nabla\zeta$','Interpreter','latex','fontsize',18)%'Rotation',0)
hold on
plot(x(45:end-45),-d(45:end-45,la),'k','LineWidth',1)
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
ylabel('Depth [m]','fontsize',14)
set(gca,'fontsize',14)
set(gca,'Layer','top')
tobj=title('(a)','fontsize',12);
tobj.Position=[x(45)+0.4,zc(1)+0.1];

subplot(4,1,2)
pcolor(X,Z,1e11*sadvsecstmonth0(45:end-44,1:168)')
shading flat
caxis([-5 5])
colormap(map)
colorbar('Ticks',-4:2:4)
% ylabel(colorbar,'\boldmath$\bf{u}\cdot\nabla\zeta$','Interpreter','latex','fontsize',14)%'Rotation',0)
ylabel(colorbar,'$\bf{u}\cdot\nabla\zeta$','Interpreter','latex','fontsize',18)%'Rotation',0)
hold on
plot(x(45:end-45),-d(45:end-45,la),'k','LineWidth',1)
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
ylabel('Depth [m]','fontsize',14)
set(gca,'fontsize',14)
set(gca,'Layer','top')
tobj=title('(b)','fontsize',12);
tobj.Position=[x(45)+0.4,zc(1)+0.1];

subplot(4,1,3)
pcolor(X,Z,1e11*bvsecwkmonth0(45:end-44,1:168)')
shading flat
caxis([-0.5 0.5])
colormap(map)
colorbar('Ticks',-0.4:0.2:0.4)
% ylabel(colorbar,'\boldmath$\beta v$','Interpreter','latex','fontsize',14,'fontweight','bold')%'Rotation',0)
ylabel(colorbar,'$\beta v$','Interpreter','latex','fontsize',20)%'Rotation',0)
hold on
plot(x(45:end-45),-d(45:end-45,la),'k','LineWidth',1)
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
ylabel('Depth [m]','fontsize',14)
set(gca,'fontsize',14)
set(gca,'Layer','top')
tobj=title('(c)','fontsize',12);
tobj.Position=[x(45)+0.4,zc(1)+0.1];
colormap(map)


subplot(4,1,4)
pcolor(X,Z,1e11*bvsecstmonth0(45:end-44,1:168)')
shading flat
caxis([-0.5 0.5])
colormap(map)
colorbar('Ticks',-0.4:0.2:0.4)
% ylabel(colorbar,'\boldmath$\beta v$','Interpreter','latex','fontsize',14,'fontweight','bold')%'Rotation',0)
ylabel(colorbar,'$\beta v$','Interpreter','latex','fontsize',20)%'Rotation',0)
hold on
plot(x(45:end-45),-d(45:end-45,la),'k','LineWidth',1)
xticks(140:5:160)
xticklabels({'140^{\circ}E','145^{\circ}E','150^{\circ}E',...
    '155^{\circ}E','160^{\circ}E'})
ylabel('Depth [m]','fontsize',14)
set(gca,'fontsize',14)
set(gca,'Layer','top')
tobj=title('(d)','fontsize',12);
tobj.Position=[x(45)+0.4,zc(1)+0.1];
colormap(map)
print -dpng fig15_sectionsense.png

















