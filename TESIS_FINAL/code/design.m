%parametros diseño Antena
clc
clear all
close all

f=12e9;
lambda=3e8/f;
t=0.3e-3;
T=3e-3;
c=lambda/4;
a=20e-2;
eps0=8.854187817e-12;
mu0=4*pi*10^-7;

bs=1.4e-3;
fprintf(sprintf('bs=%8.10g mm\n',bs*1e3))

alpha=asin((bs+t)/c)/2;
fprintf(sprintf('alpha=%8.10g°\n',rad2deg(alpha)))

d=(bs+t)/sin(alpha);
fprintf(sprintf('d=%8.10g mm\n',d*1e3))

DELTA=0.1*lambda;
fprintf(sprintf('Delta''=%8.10g mm\n',DELTA*1e3))

psi=atan(DELTA*sin(alpha)/(DELTA*cos(alpha)+20*d))
psi=acot(19*d/DELTA / sin(alpha) - cot(alpha))
fprintf(sprintf('psi''=%8.10g°\n',rad2deg(psi)))

d_=sqrt((lambda/4 *sin(2*alpha))^2 + (lambda/4 *(1+cos(2*alpha)) +DELTA/19 )^2); %correcto
fprintf(sprintf('d''=%8.10g mm\n',d_*1e3))

d_=sqrt((d*sin(alpha))^2 + (d*cos(alpha)+DELTA/19)^2); %correcto
fprintf(sprintf('d''=%8.10g mm\n',d_*1e3))

d_=(DELTA/19 + d*cos(alpha)) / cos(alpha-psi); %casi correcto
fprintf(sprintf('d''=%8.10g mm\n',d_*1e3))

d_=DELTA/19*sin(alpha) /(sin(psi));
fprintf(sprintf('d''=%8.10g mm\n',d_*1e3))
d_=DELTA/19*sin(pi-alpha)/ sin(psi) ;
fprintf(sprintf('d''=%8.10g mm\n',d_*1e3))

d_ =sqrt((20*d)^2 + DELTA^2 - 2*20*d*DELTA*cos(pi-alpha))/20;
fprintf(sprintf('d''=%8.10g mm\n',d_*1e3))

d_=(bs+t)/sin(alpha-psi);
fprintf(sprintf('d''=%8.10g mm\n',d_*1e3))

%frecuencia angular
omega=2*pi*f;
%numero de onda en el vacío
k=omega*sqrt(mu0*eps0)

%onda de corte
kx= pi/a           %mode 1 ,x dependent, TE1_
kq=0/1.4e-3  %mode n , q dependent, TE1n
            
kc = sqrt(kx^2 + kq.^2)   %Shear wave number

fc=kc/k * f
fc= 1/(2*sqrt(mu0*eps0))*sqrt((1/a)^2)
%phase constant
beta = sqrt(k^2 - kc.^2)

lambda1=2*pi/beta
lambda1= lambda/sqrt(1-(fc/f)^2)


%% iNCLINACION CON VARIACION DE C
deltac=(41.65/19)/360 * lambda;
psi = atan(deltac*sin(alpha)/(deltac*cos(alpha)+d))
d_ =deltac*sin(alpha)/sin(psi)
d_= sqrt( (deltac*sin(alpha))^2 + (deltac*cos(alpha)+d)^2)
