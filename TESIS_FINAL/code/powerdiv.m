clear all;
close all;
% clear figures;
clc;
% Cristian Duguet Sáez. 2009-09-12
% Algoritmo para el cálclo de altura de steps en un divisor de potencia squintless
% 20 elementos
% Distribución, como vector A.
% Recordar que se realiza la distribución de potencia para la mitad del
% arreglo (arreglo final tiene 40 elementos)
%freq: 12GHz;
%lambda = 25 %mm;

el=20;

%%%%%%%%%%%%-----VALORES IMPORTADOS DE DIST. TAYLOR-----------%%%%%%%%%%%%%
%n=7, SLL=30dB
A=[40.0173   39.6844   39.0521   38.1584   37.0171   35.6173   33.9545   32.0634   30.0176   27.8902   25.7070   23.4396   21.0569   18.6039   16.2431   14.2140   12.7217   11.8276   11.4199   11.2941]
%%%%%%%%%%%%--------------------------------------------------%%%%%%%%%%%%%

P=A.^2; P=P/sum(P);

Ptot=sum(P);
% P=0.995*P;

N=length(P) %nro de elementos

%% ---------------------Algoritmo de Potencias ------------------------- %%
%correccion1
% P=P-[-0.0171251292826020,-0.00143750155408756,-0.000620342805139687,0.000272357306311422,0.00135663524698268,0.00251649179597321,0.00315998908320211,0.00347218967929941,0.00320868230829512,0.00296559121089881,0.00246542071915242,0.00204891257661037,0.00154295939198814,0.000968453520189275,0.000268547843438473,-0.000429784673302779,-0.000957050569550727,-0.00130987650473384,-0.00142099666507059,-0.00137573424570286;];
% P=P-[0.0228130254494288,-0.00324636592076200,-0.00463520986089473,-0.00234560261230657,-0.00118516119946133,-0.00232419969424526,-0.00156674786827761,-0.00253317692862878,-0.00172654049412993,-0.00206365306791859,-0.00152425137161466,-0.00141143581768854,-0.00100151769169605,-0.000965085360767611,-0.000620433719929838,-0.000651302564072148,-0.000232597786400770,-0.000368813970540973,0.000138504841077983,-0.000603585002080292;];


br=1; %mm

b1min=br*(Ptot-P(1))/(P(1));
b1max=br*(Ptot-P(1))/(P(1)-P(N));
b(1)=b1min;
b0=b(1)+br

%-------------------- DETERMINACION DE b() -------------------------------%
for i=2:N
    b(i) = P(i-1)/P(i)*b(i-1);
    b(i) = b(i)-br;
end

%---------------------- NO   CRECIMIENTO DE B ----------------------------%
%  figure(1)
%  limit(1)=Ptot*br/(br+b(1));
%  for i=[2:N]
%      b(i) = min([b(i) b(i-1)]);
%      limit(i)=  (Ptot- sum(P(1:i-1))) * br/(br+b(i));
%  end
%     
% % Cómo se ve afectada la distribución por esta restricción:
%  figure(1)
%  semilogy([1:20],P,'b',[1:20],limit,'r');
%  title('Distribucion de Potencia'); legend('Distrib. de Taylor','Sin decrecimiento de steps');
% 
% clear a,i;
clear P;
%% -------------------Algoritmo de Impedancias ------------------------- %%
%normalización de potencia distribuída
P=A.^2; P=P/sum(P);
clear br;
br=0.52731082463;
br=1.4;
br=1;
% ls = 0.795492119137066; % pérdida de acoplamiento obtenida de ajuste de curvas.
% ls=0.8;
% ls=0.95
% lm=1/1.01;
%-------- energía hacia la carga----%
% P=P-[-0.0173906357852967,-0.00273850017516476,-0.00236961632448232,-0.00164550971692996,-0.000658556072828956,0.000239198145070979,0.000924636576497911,0.00173632804064307,0.00190720328118452,0.00231574929413482,0.00238332388683540,0.00246086352339619,0.00225013353477933,0.00204391152657028,0.00169456748552252,0.00138838215246226,0.00123350190570444,0.00116645903089515,0.00128748256015125,0.00148801637328809;];
% P=P-[0.00609879683329602,-0.00250623744638762,-0.00149432246403330,-0.000308645192879473,4.82221373034469e-05,1.13606703260744e-05,0.000149837555059262,-0.000275434146401257,8.78296615160687e-05,-0.000156489966915557,-4.89094084434594e-05,-0.000329671967831055,-5.84475146998692e-05,-0.000279446775222826,-2.49269620451395e-05,-0.000211952416481724,2.54083821017402e-05,-0.000260616866793779,6.85249405684377e-05,-0.000698829969121926;];
% P=P-[-0.00165759369090024,0.00135970020644888,0.00109412778991397,0.000318901432475760,-5.23018158324340e-05,-5.73402239337911e-05,-0.000119900536128412,6.14925178324538e-05,1.52670976556121e-05,-0.000139355336198549,-8.69645813735623e-05,-5.65810931769523e-05,-6.82805611795512e-05,-0.000158728208503978,-9.22978880109511e-05,-0.000121555165690965,-6.45725889897083e-05,-0.000121884415658434,-2.41293458562945e-05,-0.000223748282131222;];
% P=0.995*P; %Para dejar un remanente a la carga
%----------------------------------%
cumP=cumsum(P);
%Vector de acoplamientos:
for i=1:N
    if (i==1) S2(i) = P(i); 
    else      S2(i) = P(i) / (1 - cumP(i-1));
    end
end
%  S2 =S2-correccion-correccion2;

% b2=zeros(size(P)); %se inicializa vector de alturas  para el nuevo algoritmo

for i=[1:1:N-1]
    b2(i) = br* sqrt((1-S2(i))/(S2(i)*S2(i+1)));
end
b2(N) = br* sqrt((1-S2(N))/(S2(N))) ;    %Se considera que S2(N+1)=1 = full absorción por la carga

b02=br/S2(1);
 
%RESULTADOS ARNDT
ARNDT=1e3*[0.0165487878173164,0.0153978770396513,0.0143305129211123,0.0133800579892415,0.0124991858255087,0.0116765179805423,0.0109275072336703,0.0102328200886909,0.00959265652681797,0.00896890601397655,0.00834292308054726,0.00771249574176216,0.00708556643697299,0.00648006612970983,0.00589966905183740,0.00531085259444016,0.00464091737627374,0.00379902585716515,0.00272613426905429,0.00142638452377978,1.88451286189225e-05;];
  
%% Gráfico
figure(1);
step1 = -[b0 b0  b b(N)];
step2 = -[b02 b02 b2 b2(N)];
step3 = -[ARNDT(1) ARNDT ARNDT(length(ARNDT))];
stairs([-1:21],[step1 ; step2]'); grid on;
legend('Algoritmo de Potencias','Algoritmo de Impedancias');hold on;
plot([-1:21],zeros(length(step1)),'r');
ax=axis; ax(1:2) =[0 21]; ax(4)=1;
axis(ax);
title('Comparación de algoritmos');
ylabel('Altura de step (mm)');
xlabel('step');
%Exportar Gráfico
set(gcf, 'PaperUnits', 'centimeters');
myfiguresize=[0 0 17 9]; 
set(gcf, 'PaperPosition', myfiguresize);
print -dpng -f1 -r300 desarrollo_alturasteps
%% --------------------------------------------------------------------- %%
%% --------------------------------------------------------------------- %%


fout = fopen('completedist_results/steps_cambiar_12GHz.vbs','wt');
% fout = fopen('completedist_results/steps_cambiar_12GHz_inclinado.vbs','wt');

fprintf(fout,''' ---------------------------------------------- \n');
fprintf(fout,''' Script Recorded by Cristian Duguet S. \n');
fprintf(fout,''' using MATLAB. Run it in Ansoft HFSS \n');
fprintf(fout,''' ---------------------------------------------- \n');
fprintf(fout,'Dim oAnsoftApp \n');
fprintf(fout,'Dim oDesktop \n');
fprintf(fout,'Dim oProject \n');
fprintf(fout,'Dim oDesign \n');
fprintf(fout,'Dim oEditor \n');
fprintf(fout,'Dim oModule \n');
fprintf(fout,'Set oAnsoftApp = CreateObject("AnsoftHfss.HfssScriptInterface") \n');
fprintf(fout,'Set oDesktop = oAnsoftApp.GetAppDesktop() \n');
fprintf(fout,'oDesktop.RestoreWindow \n');
fprintf(fout,'Set oProject = oDesktop.SetActiveProject("TEM_divider2") \n');
%fprintf(fout,'Set oDesign = oProject.SetActiveDesign("steps 12GHz") \n');
fprintf(fout,'Set oDesign = oProject.SetActiveDesign("steps square with absorber") \n');
fprintf(fout,'oDesign.ChangeProperty Array("NAME:AllTabs", Array("NAME:LocalVariableTab", Array("NAME:PropServers",  _ \n');
fprintf(fout,'  "LocalVariables"), Array("NAME:ChangedProps", _ \n');
fprintf(fout,' Array("NAME:b1", "Value:=", "%3.6gmm"), _ \n',b2(1));
fprintf(fout,' Array("NAME:b2", "Value:=", "%3.6gmm"), _ \n',b2(2));
fprintf(fout,' Array("NAME:b3", "Value:=", "%3.6gmm"), _ \n',b2(3));
fprintf(fout,' Array("NAME:b4", "Value:=", "%3.6gmm"), _ \n',b2(4));
fprintf(fout,' Array("NAME:b5", "Value:=", "%3.6gmm"), _ \n',b2(5));
fprintf(fout,' Array("NAME:b6", "Value:=", "%3.6gmm"), _ \n',b2(6));
fprintf(fout,' Array("NAME:b7", "Value:=", "%3.6gmm"), _ \n',b2(7));
fprintf(fout,' Array("NAME:b8", "Value:=", "%3.6gmm"), _ \n',b2(8));
fprintf(fout,' Array("NAME:b9", "Value:=", "%3.6gmm"), _ \n',b2(9));
fprintf(fout,'Array("NAME:b10", "Value:=", "%3.6gmm"), _ \n',b2(10));
fprintf(fout,'Array("NAME:b11", "Value:=", "%3.6gmm"), _ \n',b2(11));
fprintf(fout,'Array("NAME:b12", "Value:=", "%3.6gmm"), _ \n',b2(12));
fprintf(fout,'Array("NAME:b13", "Value:=", "%3.6gmm"), _ \n',b2(13));
fprintf(fout,'Array("NAME:b14", "Value:=", "%3.6gmm"), _ \n',b2(14));
fprintf(fout,'Array("NAME:b15", "Value:=", "%3.6gmm"), _ \n',b2(15));
fprintf(fout,'Array("NAME:b16", "Value:=", "%3.6gmm"), _ \n',b2(16));
fprintf(fout,'Array("NAME:b17", "Value:=", "%3.6gmm"), _ \n',b2(17));
fprintf(fout,'Array("NAME:b18", "Value:=", "%3.6gmm"), _ \n',b2(18));
fprintf(fout,'Array("NAME:b19", "Value:=", "%3.6gmm"), _ \n',b2(19));
fprintf(fout,'Array("NAME:b20", "Value:=", "%3.6gmm"), _ \n',b2(20));
fprintf(fout,'  Array("NAME:B", "Value:=", "%3.6gmm")))) ',b02);
fclose(fout);