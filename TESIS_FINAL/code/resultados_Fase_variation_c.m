%Cristian Duguet
%2009-10-14
%Resultados de simulación en divisor de potencia.
%Uniformidad de Fase en el divisor.

clear all;
close all;
clc

%---------------------IMPORTAR DATOS--------------------------------------%          
newData1 = importdata('completedist_results/phase_variation_c.csv');

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

%--------------------extraccion de variables------------------------------%
%frecuencia
freq = data(:,1);

%reflexión:
phi_11(:,i)=deg2rad(data(:,1+i));

for i=1:20    
    %parámetros S
    phi(:,i,:)=deg2rad(data(:,7+5*(i-1):11+5*(i-1)));
end


%% -----------------------------------------------------------------------%
%hacer continuos curvas que sobrepasan el límite [-pi,+pi]
phi=rad2deg(unwrap(phi,[],2));

%-------------------------------------------------------------------------%
%------------------UNIFORMIDAD DE FASE -----------------------------------%
for i=1:5
figure(i)
plot([1:20],phi(:,:,i));
xlabel('elemento'); 
ylabel('ángulo (°)');
end

%-------------------------------------------------------------------------%
%VARIANZA

vars = var(phi');
std=sqrt(vars);
% ERRORES=sqrt(sum((b-S).^2,2)/size(S,2));
% ERROR=sqrt(sum(sum((b-S).^2))/prod(size(S)))

std = phi(:,20)-phi(:,1);

% 
% figure(2)
% set(gcf,'Position',[20 100 700 600])
%  subplot(4,3,1);plot([1:20],phi(1,:));  xlabel(sprintf('f=10.2GHz,       \\sigma=%8.6g',std(1)));
%  subplot(4,3,2);plot([1:20],phi(2,:));  xlabel(sprintf('f=10.5GHz,       \\sigma=%8.6g',std(2)));title('Resultados Distribuidor Inclinado');
%  subplot(4,3,3);plot([1:20],phi(3,:));  xlabel(sprintf('f=10.8GHz,       \\sigma=%8.6g',std(3)));
%  subplot(4,3,4);plot([1:20],phi(4,:));  xlabel(sprintf('f=11.1GHz,       \\sigma=%8.6g',std(4)));
%  subplot(4,3,5);plot([1:20],phi(5,:));  xlabel(sprintf('f=11.4GHz,       \\sigma=%8.6g',std(5)));
%  subplot(4,3,6);plot([1:20],phi(6,:));  xlabel(sprintf('f=11.7GHz,       \\sigma=%8.6g',std(6)));
%  subplot(4,3,7);plot([1:20],phi(7,:));  xlabel(sprintf('f=12GHz,         \\sigma=%8.6g',std(7)));
%  subplot(4,3,8);plot([1:20],phi(8,:));  xlabel(sprintf('f=12.3GHz,       \\sigma=%8.6g',std(8)));
%  subplot(4,3,9);plot([1:20],phi(9,:));  xlabel(sprintf('f=12.6GHz,       \\sigma=%8.6g',std(9)));
%  subplot(4,3,10);plot([1:20],phi(10,:));xlabel(sprintf('f=12.9GHz,       \\sigma=%8.6g',std(10)));
%  subplot(4,3,11);plot([1:20],phi(11,:));xlabel(sprintf('f=13.2GHz,       \\sigma=%8.6g',std(11)));
%  subplot(4,3,12);plot([1:20],phi(12,:));xlabel(sprintf('f=13.5GHz,       \\sigma=%8.6g',std(12)));
% %  subplot(4,3,13);plot([1:20],phi(13,:));xlabel(sprintf('f=13.8GHz,       \\sigma=%8.6g',std(13)));



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
%-----------------COMPARACION ENTRE 11.1GHZ-13.5GHZ-------------------------------------------%
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
% 
%-------------------------------------------------------------------------%