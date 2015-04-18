clear all
clc 
close all

% Parámetros de diseño del reflector pillbox parabólico
%Cristian Duguet 

f=12e9;
lambda=3e8/f;
%ancho de la guía resultante
a=20e-2; 
%se tomará el mismo ancho para diseñar el tamaño del reflector??

%Relacion f/D con ángulo máximo : tan(beta/2)= a/(2*f);
beta=deg2rad(90); %degrees
a=20e-2; %cm

f = a/(2*tan(beta/2))

%altura de la forma parabólica: 
h= f/(16*(f/(2*a))^2)



%Feedhorn
%H plane horn

W=0.0256;
BW=50.6 * lambda/W;

a0=1.905e-2;
L=1.5e-2;
LH = W/((W-a0)/L);
deltaphi=rad2deg(W^2/(8*lambda*LH));


obstruccion= W/2*cos(deg2rad(50)); %lo que se le quita a la apertura
D_2=20+obstruccion;
f=D_2/2;


% Archivo de resultados de Pillbox
%Distribucion de apertura 
PillBoxData = importdata('pillbox_distribution.csv');

distrib=PillBoxData.data(:,2);

dist = linspace(0,200,length(distrib));

figure(1)
plot(dist,distrib); title('Distribución de Amplitud PillBox');
xlabel('length (mm)'); ylabel('E field complex magnitude');
grid on; set(1, 'PaperUnits', 'centimeters');
set(1, 'PaperPosition', [0 0 17 9]);
print -dpng -f1 -r300 desarrollo_pillbox_distribution


pattern= fftshift(fft(dist,length(dist)));
% pattern=fftshift(pattern)
figure(2)
plot(dist,pattern)


