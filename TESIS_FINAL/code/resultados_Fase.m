%Cristian Duguet
%2009-10-14
%Resultados de simulación en divisor de potencia.
%Uniformidad de Fase en el divisor.

clear all;
close all;
clc

%---------------------IMPORTAR DATOS--------------------------------------%          
newData1 = importdata('completedist_results/phase_br4.20_variation_c.csv');

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

%--------------------extraccion de variables------------------------------%
%frecuencia
freq = data(:,1);
%reflexión:
phi_11 =deg2rad(data(:,2));
%parámetros S
phi=deg2rad(data(:,3:22));


%% -----------------------------------------------------------------------%
%hacer continuos curvas que sobrepasan el límite [-pi,+pi]
phi=rad2deg(unwrap(phi,[],2));

%-------------------------------------------------------------------------%
%------------------UNIFORMIDAD DE FASE -----------------------------------%
figure(1)
plot([1:20],phi);
xlabel('elemento'); 
ylabel('ángulo (°)');
% % daspect([100 50 1]);
% legend('Distribución de Taylor-7','Simulación');
% ax=axis

%-------------------------------------------------------------------------%
%Error de ángulo
Recta =@(param,x)param(1)*x+param(2);
xdata=1:20;
for i=1:13
ydata=phi(i,:);
%AJUSTE CON MINIMOS CUADRADOS
[x,resnorm] = lsqcurvefit(Recta,[0 phi(1,1)], xdata, ydata);
a(i)=x(1); b(i)=x(2);
end

%No Uniformidad de fase
for i=1:13
approx(i,:)=a(i)*[1:20]+b(i);
RIZADO(i)=std(phi(i,:)-approx(i));
end

std = a*19;

PROMEDIODESFASE=mean(std)

figure(2)
set(gcf,'Position',[20 100 700 600])
 subplot(4,3,1);plot([1:20],phi(1,:),[1:20],approx(1,:),'--');  xlabel(sprintf( 'f=10.2GHz,     \\Delta \\phi =%8.6g',std(1)));
 subplot(4,3,2);plot([1:20],phi(2,:),[1:20],approx(2,:),'--');  xlabel(sprintf( 'f=10.5GHz,     \\Delta \\phi =%8.6g',std(2)));title('Uniformidad de Fase');
 subplot(4,3,3);plot([1:20],phi(3,:),[1:20],approx(3,:),'--');  xlabel(sprintf( 'f=10.8GHz,     \\Delta \\phi =%8.6g',std(3)));
 subplot(4,3,4);plot([1:20],phi(4,:),[1:20],approx(4,:),'--');  xlabel(sprintf( 'f=11.1GHz,     \\Delta \\phi =%8.6g',std(4)));
 subplot(4,3,5);plot([1:20],phi(5,:),[1:20],approx(5,:),'--');  xlabel(sprintf( 'f=11.4GHz,     \\Delta \\phi =%8.6g',std(5)));
 subplot(4,3,6);plot([1:20],phi(6,:),[1:20],approx(6,:),'--');  xlabel(sprintf( 'f=11.7GHz,     \\Delta \\phi =%8.6g',std(6)));
 subplot(4,3,7);plot([1:20],phi(7,:),[1:20],approx(7,:),'--');  xlabel(sprintf( 'f=12GHz,       \\Delta \\phi =%8.6g',std(7)));
 subplot(4,3,8);plot([1:20],phi(8,:),[1:20],approx(8,:),'--');  xlabel(sprintf( 'f=12.3GHz,     \\Delta \\phi =%8.6g',std(8)));
 subplot(4,3,9);plot([1:20],phi(9,:),[1:20],approx(9,:),'--');  xlabel(sprintf( 'f=12.6GHz,     \\Delta \\phi =%8.6g',std(9)));
 subplot(4,3,10);plot([1:20],phi(10,:),[1:20],approx(10,:),'--');xlabel(sprintf('f=12.9GHz,     \\Delta \\phi =%8.6g',std(10)));
 subplot(4,3,11);plot([1:20],phi(11,:),[1:20],approx(11,:),'--');xlabel(sprintf('f=13.2GHz,     \\Delta \\phi =%8.6g',std(11)));
 subplot(4,3,12);plot([1:20],phi(12,:),[1:20],approx(12,:),'--');xlabel(sprintf('f=13.5GHz,     \\Delta \\phi =%8.6g',std(12)));
%  subplot(4,3,13);plot([1:20],phi(13,:),[1:20],approx(13,:),'--');xlabel(sprintf('f=13.8GHz,     \\Delta \\phi =%8.6g',std(13)));

set(2, 'PaperUnits', 'centimeters');
myfiguresize=[0 0 17 20]; 
set(2, 'PaperPosition', myfiguresize);
print -dpng -f2 -r300 desarrollo_fase_frecuencias

% 
% figure(4)
% plot(10.2:0.3:13.8,mag2db(S_11)); title('Pérdidas de Retorno')
% xlabel('Frecuencia (GHz)'); ylabel('dB');




% export settings
% figure(1)
% set(1, 'PaperUnits', 'centimeters');
% myfiguresize=[0 0 10 7]; 
% set(1, 'PaperPosition', myfiguresize);
% ax=axis; ax(3:4)=[0 0.6]; axis(ax);
% print -dpng -f1 -r300 desarrollo_completedist_results
% 
% set(2, 'PaperUnits', 'centimeters');
% myfiguresize=[0 0 17 9]; 
% set(2, 'PaperPosition', myfiguresize);
% print -dpng -f2 -r300 desarrollo_distribuciones_algoritmos

% set(3, 'PaperUnits', 'inches');
% myfiguresize=[0 0 7 9]; 
% set(3, 'PaperPosition', myfiguresize);
% print -dpng -f3 -r300 desarrollo_completedist_frequencies


% set(4, 'PaperUnits', 'centimeters');
% myfiguresize=[0 0 10 7]; 
% set(4, 'PaperPosition', myfiguresize);
% print -dpng -f4 -r300 desarrollo_completedist_reflection
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-----------------COMPARACION ENTRE 10.8GHZ-12.9GHZ-------------------------------------------%
% ERROR CUADRÁTICO MEDIO
% [a,b]= ndgrid(ones(9,1),Taylor);
% ERROR2=sqrt(sum(sum((b-S(4:12,:)).^2))/(9*20))
% 
% figure(5)
% plot([1:20],Taylor,'r-x',[1:20],S(4:12,:));
% % hold on; for i=1:20
% title('Distribución de Apertura');
% xlabel('elemento'); 
% ylabel('Voltaje');
% % daspect([100 50 1]);
% legend('Distribución de Taylor-7','Simulación');
% 
% set(5, 'PaperUnits', 'centimeters');
% myfiguresize=[0 0 10 7]; 
% set(5, 'PaperPosition', myfiguresize);
% % ax=axis; ax(3:4)=[0 0.6]; axis(ax);
% print -dpng -f5 -r300 desarrollo_RESULTS_2
% 
figure(6)
plot([1:20],phi(3:10,:)); title('Distribución de Fase en el Divisor');grid on;
xlabel('elemento'); ylabel('ángulo (°)'); legend('10.8GHz','11.1GHz','11.4GHz','11.7GHz','12.0GHz','12.3GHz','12.6GHz','12.9GHz');

set(6, 'PaperUnits', 'centimeters');
myfiguresize=[0 0 17 9]; 
set(6, 'PaperPosition', myfiguresize);
print -dpng -f6 -r300 evaluacion_fase

% DESFASE=std(3:10)
PROMEDIODESFASE=mean(std(3:10)); %Desfase en el rango de interés
%-------------------------------------------------------------------------%
