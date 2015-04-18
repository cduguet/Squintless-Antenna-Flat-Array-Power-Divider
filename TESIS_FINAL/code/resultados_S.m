%Cristian Duguet
%2009-10-14
%Resultados de simulación en divisor de potencia.
%parametros S en HFSS.

clear all;
close all;
clc

%---------------------IMPORTAR DATOS--------------------------------------%          
%Datos de simulació importados  en /completedist_results/steps_generic_091014.csv (para div cuadrado)
newData1 = importdata('completedist_results/p100_br4.20_variation_c.csv');
% newData1 = importdata('completedist_results/square/square_p100_it1_abs_graphite.csv');


% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

% open('completedist_results/steps_generic_091015.mat');
%Datos de simulación importados  en S_sloping.mat (para div inclinado)
%open('S_sloping.mat');
            %--------------------------------------%
%--------------------extraccion de variables------------------------------%
%frecuencia
freq = data(:,1);
%reflexión:
S_11 =data(:,2);
%parámetros S
S=data(:,3:22);


%% Resultados  10GHz %%%%%%%%%%5
% Algoritmo Chtcherbakov, 08-10-2009, 10GHz
% S=[0.355736082416702,0.359235969914973,0.346713548416398,0.331682538866723,0.314195893624868,0.293995897456904,0.271915150917994,0.248126157009361,0.223845936556941,0.198949265984766,0.175682847440354,0.151432209553191,0.129226600158105,0.106411689398144,0.0878656780824335,0.0704652914130182,0.0603719869983563,0.0493849759435185,0.0466930101098820,0.0256879291654442;]
% % % % % 
% % % % % %Parámetros obtenidos con método Duguet
% S=[0.316696414684651,0.334201674194120,0.328914073701039,0.320538634942077,0.308210726854087,0.295283842698755,0.277614503252343,0.259493400577042,0.237904362513122,0.217914351567008,0.194498558584651,0.174822645140965,0.150794957722852,0.130511013701462,0.107661828624655,0.0908483282383447,0.0756348560873854,0.0675958943410470,0.0610295189271969,0.0635394195161871;];
% % Parámetros obtenidos con el algoritmo chtcherbakov simple, 2010-01-28, 10GHz
% S_chtcherbakov=[0.328573932603124,0.342126944339365,0.332319391698228,0.321925357473024,0.307268651885562,0.292619672827679,0.273138932885255,0.254173226630989,0.232412890590959,0.212270824346091,0.189898336506527,0.169469527616838,0.147334409242359,0.126522342856417,0.107411545016724,0.0901877351177687,0.0799967483478288,0.0700392986055468,0.0691962251701893,0.0643653846601728;];
% S_ref_chtcherbakov = 0.013735056805785;
% % Parámetros obtenidos con el algoritmo duguet simple, 2010-01-28,% 10GHz
% S_duguet=[0.249888569514298,0.321266569925012,0.312108597419067,0.304709000877137,0.296005598482215,0.290592588053904,0.276497482980507,0.266160790702200,0.247025445548784,0.233611717922029,0.212460578313455,0.198156117080089,0.174172335655313,0.157288017617940,0.132894464488890,0.116422262930646,0.100598935284635,0.0937296002992690,0.0880044205426786,0.102535588929884;]
% S_ref_duguet=[0.0356747468231971;]

%% Resultados 12GHz
%Chtcherbakov = Algoritmo de Impedancias Adaptadas
%Duguet = Algoritmo de Razón de Potencias

%Potencia 100%
S_impedancias = [0.291082427588308,0.312552613488811,0.308031241089217,0.302005062474192,0.294487217461738,0.284845203348659,0.272846308561320,0.259413555847595,0.243653836501606,0.227859096580931,0.211010534912138,0.193640960276058,0.174714408202275,0.155289698678899,0.136086100560688,0.119466518184261,0.107489893230023,0.100435660319193,0.0979996779611993,0.0980930890723278;];
S_impedancias2 = [0.328965670350715,0.312923950709215,0.309448766918552,0.304210325735487,0.295684797476112,0.284444990105282,0.271422749735722,0.255506615903481,0.239891263747313,0.222368002719878,0.205166304613179,0.186296231860461,0.167977806224489,0.147619552275402,0.129614552879860,0.112569598018854,0.101715208415453,0.0930604425417341,0.0915695323834077,0.0862288106223298;];

S_potencias = [0.291538138419997,0.314626977259998,0.310857715050320,0.305163767153445,0.297889094409595,0.288814964158094,0.276912369897488,0.262737996102866,0.246310111585095,0.229280635491622,0.211204977868003,0.192574324742028,0.172678748810269,0.151787458340280,0.130741757383514,0.111597859035491,0.0967653071681508,0.0872409670288181,0.0830388924256395,0.0822101180193818;];
S_potencias2 = [0.353458117583156,0.311739106390406,0.304331483666332,0.300843754903877,0.293591750772607,0.280309457620364,0.268241912554253,0.251049572770442,0.236079326216077,0.218037532394367,0.201538757030643,0.183369905261330,0.165146823181664,0.145278675752656,0.127296604669415,0.110600923371082,0.100438924001089,0.0924772883594464,0.0919508518798249,0.0867793336479559;];
S_potencias3 = [0.247898187276580,0.327261126636867,0.322600558440901,0.319748294279611,0.293543299903729,0.295951992845697,0.271223575655706,0.265122347420211,0.241750067889139,0.229599757728416,0.208199450199620,0.192643633222848,0.170322580527900,0.152516198917418,0.131466355376442,0.115502916957178,0.104553988770270,0.0919371535773555,0.101889773745852,0.0389441034652087;];

%% -------------------------------------------------------------------------%
%Distribución de Taylor Teórica,  SLL=30dB, n_=7
Taylor=[40.0173   39.6844   39.0521   38.1584   37.0171   35.6173   33.9545   32.0634   30.0176   27.8902   25.7070   23.4396   21.0569   18.6039   16.2431   14.2140   12.7217   11.8276   11.4199   11.2941];
Taylor=Taylor/norm(Taylor); %Voltaje V
%-------------------------------------------------------------------------%
%------------------COMPARACIÓN DE DISTRIBUCIONES--------------------------%
% figure(1)
% plot([1:20],Taylor,'r-x',[1:20],S);
% % hold on; for i=1:20
% title('Distribución de Apertura');
% xlabel('elemento'); 
% ylabel('Voltaje');
% % daspect([100 50 1]);
% legend('Distribución de Taylor-7','Simulación');
% ax=axis

%-------------------------------------------------------------------------%
%ERROR

% ERROR CUADRÁTICO MEDIO
[a,b]= ndgrid(ones(13,1),Taylor);
ERRORES=sqrt(sum((b-S).^2,2)/size(S,2));
ERROR=sqrt(sum(sum((b-S).^2))/prod(size(S)))

%Factor de Corrección
Delta=S(7,:).^2-Taylor.^2;

figure(2)
plot([1:20],Taylor,'r-x',[1:20],S_impedancias,'-.',[1:20],S_potencias,'m');
% hold on; for i=1:20
title('Comparación de distribuciones de amplitud.');
xlabel('elemento'); 
ylabel('Voltaje'); grid on;
% daspect([100 1 1]);
legend('Distribución de Taylor-7','Algoritmo de Impedancias','Algoritmo de Razón de Potencias');


% figure(2)
% [a,b] = ndgrid(8.5:.3:12.4,1:20);
% surface(a,b,Scuad);xlabel('f');ylabel('elemento');

% figure(3)
% subplot(4,3,1);plot([1:20],Taylor,'r-',[1:20],S(1,:));  xlabel(sprintf('f=10.2GHz,       \\sigma=%8.6g',ERRORES(1)));ax=axis(); ax(3:4)=[0.05 0.35];axis(ax);
% subplot(4,3,2);plot([1:20],Taylor,'r-',[1:20],S(2,:));  xlabel(sprintf('f=10.5GHz,       \\sigma=%8.6g',ERRORES(2)));title('Resultados Distribuidor Inclinado');axis(ax);
% subplot(4,3,3);plot([1:20],Taylor,'r-',[1:20],S(3,:));  xlabel(sprintf('f=10.8GHz,       \\sigma=%8.6g',ERRORES(3)));axis(ax);
% subplot(4,3,4);plot([1:20],Taylor,'r-',[1:20],S(4,:));  xlabel(sprintf('f=11.1GHz,       \\sigma=%8.6g',ERRORES(4)));axis(ax);
% subplot(4,3,5);plot([1:20],Taylor,'r-',[1:20],S(5,:));  xlabel(sprintf('f=11.4GHz,       \\sigma=%8.6g',ERRORES(5)));axis(ax);
% subplot(4,3,6);plot([1:20],Taylor,'r-',[1:20],S(6,:));  xlabel(sprintf('f=11.7GHz,       \\sigma=%8.6g',ERRORES(6)));axis(ax);
% subplot(4,3,7);plot([1:20],Taylor,'r-',[1:20],S(7,:));  xlabel(sprintf('f=12GHz,         \\sigma=%8.6g',ERRORES(7)));axis(ax);
% subplot(4,3,8);plot([1:20],Taylor,'r-',[1:20],S(8,:));  xlabel(sprintf('f=12.3GHz,       \\sigma=%8.6g',ERRORES(8)));axis(ax);
% subplot(4,3,9);plot([1:20],Taylor,'r-',[1:20],S(9,:));  xlabel(sprintf('f=12.6GHz,       \\sigma=%8.6g',ERRORES(9)));axis(ax);
% subplot(4,3,10);plot([1:20],Taylor,'r-',[1:20],S(10,:));xlabel(sprintf('f=12.9GHz,       \\sigma=%8.6g',ERRORES(10)));axis(ax);
% subplot(4,3,11);plot([1:20],Taylor,'r-',[1:20],S(11,:));xlabel(sprintf('f=13.2GHz,       \\sigma=%8.6g',ERRORES(11)));axis(ax);
% subplot(4,3,12);plot([1:20],Taylor,'r-',[1:20],S(12,:));xlabel(sprintf('f=13.5GHz,       \\sigma=%8.6g',ERRORES(12)));axis(ax);
% subplot(4,3,13);plot([1:20],Taylor,'r-',[1:20],S(13,:));xlabel(sprintf('f=13.8GHz,       \\sigma=%8.6g',ERRORES(13)));axis(ax);



figure(4) 
plot(10.2:0.3:13.8,mag2db(S_11)); title('Pérdidas de Retorno')
xlabel('Frecuencia (GHz)'); ylabel('dB'); grid on;


% S_11= [-16.126439 -17.141122 -19.028338 -21.950059 -26.561503 -35.984942 -37.903613 -28.892389 -25.347027 -23.567148 -22.800859 -22.811909 -23.558145];
% plot(10.2:0.3:13.8,(S_11)); title('Pérdidas de Retorno')

% 
% set(4, 'PaperUnits', 'centimeters');
% myfiguresize=[0 0 17 9]; 
% set(4, 'PaperPosition', myfiguresize);
% print -dpng -f4 -r300 desarrollo_completedist_reflection
% 
% 

% export settings
% figure(1)
% set(1, 'PaperUnits', 'centimeters');
% myfiguresize=[0 0 10 7]; 
% set(1, 'PaperPosition', myfiguresize);
% ax=axis; ax(3:4)=[0 0.6]; axis(ax);
% print -dpng -f1 -r300 desarrollo_completedist_results
% 
set(2, 'PaperUnits', 'centimeters');
myfiguresize=[0 0 17 9]; 
set(2, 'PaperPosition', myfiguresize);
print -dpng -f2 -r300 desarrollo_distribuciones_algoritmos

% set(3, 'PaperUnits', 'inches');
% myfiguresize=[0 0 7 9]; 
% set(3, 'PaperPosition', myfiguresize);
% print -dpng -f3 -r300 desarrollo_completedist_frequencies


set(4, 'PaperUnits', 'centimeters');
myfiguresize=[0 0 17 9]; 
set(4, 'PaperPosition', myfiguresize);
print -dpng -f4 -r300 evaluacion_reflection
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-----------------COMPARACION ENTRE 10.8GHZ-13.5GHZ-------------------------------------------%
% ERROR CUADRÁTICO MEDIO
[a,b]= ndgrid(ones(9,1),Taylor);
ERROR2=sqrt(sum(sum((b-S(4:12,:)).^2))/(9*20))

figure(5)
plot([1:20],Taylor,'r-x',[1:20],S(4:12,:));
% hold on; for i=1:20
title('Distribución de Apertura');
xlabel('elemento'); 
ylabel('Voltaje');
% daspect([100 50 1]);
legend('Distribución de Taylor-7','Simulación');

set(5, 'PaperUnits', 'centimeters');
myfiguresize=[0 0 10 7]; 
set(5, 'PaperPosition', myfiguresize);
% ax=axis; ax(3:4)=[0 0.6]; axis(ax);
print -dpng -f5 -r300 desarrollo_RESULTS_2


[a,b]= ndgrid(ones(8,1),Taylor);
ERROR2=sqrt(sum(sum((b-S(3:10,:)).^2))/(8*20))

figure(6)
plot([1:20],Taylor,'r-x',[1:20],S(3:10,:));
% hold on; for i=1:20
title('Distribución de Apertura, frecuencia 10.8GHz-12.9GHz');
xlabel('elemento'); 
ylabel('Voltaje');
% daspect([100 50 1]);
legend('Distribución de Taylor-7','Simulación');

set(6, 'PaperUnits', 'centimeters');
myfiguresize=[0 0 17 9]; 
set(6, 'PaperPosition', myfiguresize);
% ax=axis; ax(3:4)=[0 0.6]; axis(ax);
print -dpng -f6 -r300 evaluacion_RESULTS


%-------------------------------------------------------------------------%