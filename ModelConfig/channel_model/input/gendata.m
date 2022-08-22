clear

%--------
nx = 192;
ny = 192;
nz =  40;

dx = 10.42e+3;
dy = 10.42e+3;

dz = 10 * 1.082.^(1:nz)'; 
%--------

%--------
Lx = nx*dx;
Ly = ny*dy;
Hz = sum(dz);

x  = (dx/2:dx:Lx-dx/2);
y  = (dy/2:dy:Ly-dy/2);
z  = - cumsum(dz) + dz/2;

[X,Y,Z] = ndgrid(x,y,z);
%--------

%--------
hbot = -Hz * ones(nx,ny);

hbot = hbot + 1000*( exp(-(X(:,:,1)-Lx/2).^2/(100e+3)^2) );

hbot(:,1) = 0;
%--------

%--------
taux = 0.1 * sin(pi*Y(:,:,1)/Ly);

trlx = 8 * Y(:,:,1)/Ly;

hflx = -10 * cos(3*pi*Y(:,:,1)/Ly);
hflx(Y(:,:,1)>5*Ly/6) = 0;
%--------

%--------
rbcs_temp = 8 * (exp(Z/1e+3) - exp(-Hz/1e+3)) / (1 - exp(-Hz/1e+3));

rbcs_mask = (Y-0.95*Ly)/(0.05*Ly);
rbcs_mask(Y<0.95*Ly) = 0;
%--------

%--------
T = rbcs_temp;
%--------

%--------
fid = fopen('T.init','w','b'); 
fwrite(fid,T,'real*8'); 
fclose(fid);

fid = fopen('taux.init','w','b'); 
fwrite(fid,taux,'real*8'); 
fclose(fid);

fid = fopen('trlx.init','w','b'); 
fwrite(fid,trlx,'real*8'); 
fclose(fid);

fid = fopen('hflx.init','w','b'); 
fwrite(fid,hflx,'real*8'); 
fclose(fid);

fid = fopen('rbcs_temp.init','w','b'); 
fwrite(fid,rbcs_temp,'real*8'); 
fclose(fid);

fid = fopen('rbcs_mask.init','w','b'); 
fwrite(fid,rbcs_mask,'real*8'); 
fclose(fid);

fid = fopen('topog.init','w','b'); 
fwrite(fid,hbot,'real*8'); 
fclose(fid);

fid=fopen('delZ.init','w','b'); 
fwrite(fid,dz,'real*8'); 
fclose(fid);
%--------

