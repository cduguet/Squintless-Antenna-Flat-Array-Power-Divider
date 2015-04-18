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
%freq: 10GHz;
%lambda = 30 %mm;

el=20;

%---dist coseno cuadrado--------%
%   x=0:(pi/2)/(el-1):pi/2;
% A=0.1+0.9.*(cos(x)).^2;
%  A=-x+pi/2;
%-------------------------------%

%%%%%%%%%%%%-----VALORES IMPORTADOS DE DIST. TAYLOR-----------%%%%%%%%%%%%%
%n=20, SLL=30dB
% A=[38.8482   38.5671   38.0105   37.1863   36.1113   34.7999   33.2792   31.5677   29.7022  27.7007   25.6106   23.4418   21.2570   19.0430   16.8973   14.7367   12.7930   10.6124    10.4550   19.3801] ;

%n=7, SLL=30dB
A=[40.0173   39.6844   39.0521   38.1584   37.0171   35.6173   33.9545   32.0634   30.0176   27.8902   25.7070   23.4396   21.0569   18.6039   16.2431   14.2140   12.7217   11.8276   11.4199   11.2941]
%%%%%%%%%%%%--------------------------------------------------%%%%%%%%%%%%%

P=A.^2;
P=P/sum(P);
Ptot=sum(P);

b=zeros(size(P)); %se inicializa vector de alturas de guía
N=size(P,2) %nro de elementos

clear x;


%% ---------------------Algoritmo de Duguet---------------------------- %%
%Analisis de P: para evitar decaimientos rápidos en distribución de
%potencia (más que exponencial), de tal forma que los steps siempre vayan
%aumentando. 
% 
% b0=13; %mm 
 br=1.4; %mm
% 
% % l=0.795492119137066;
% l_coef = [-0.011332382650739   0.984082130023733];
% l = @(b1)(l_coef(2) + l_coef(1)*b1);
% 
% alpha_coef = [  0.000944764920949  -0.022175172579838   0.072442771504602   0.938524216332454];
% alpha_coef = [  0.002165005287783  -0.004832942807999  -0.009174108407385   1.027863865647605 ];
% alpha = @(b1,b)(alpha_coef(1)*(b-b1).^3+alpha_coef(2)*(b-b1).^2+alpha_coef(3)*(b-b1)+alpha_coef(4));
% 
% 
% b1min = @(b1)(br*l(b1)*(Ptot-alpha(b1,b0)*P(1))/(alpha(b1 ,b0)*P(1))  -   b1);
% b1max= @(b1)(br*l(b1)*(Ptot-alpha(b1,b0)*P(1))/(alpha(b1 ,b0)*P(1)-P(N))  -   b1);
% 
% b(1)=fzero(b1max,9)
% maximo=fzero(b1max,9)

% b1min=br *l* (Ptot-P(1))/(*P(1));
% b1max=br *l* (Ptot-P(1))/(P(1)-P(N));
% b(1)=b1min+0.05;

%-------------------- DETERMINACION DE b() -------------------------------%
% for i=2:N
%     aux = @(bi)( P(i-1)/P(i) .* b(i-1)/alpha(bi,b(i-1))- l(bi)*br - bi );
%     b(i) = fzero(aux, b(i-1));
%     alpha(b(i),b(i-1))
%     if i==5
%         figure(1)
%         x = 1:.2:10
%         for k=1:size(x,2)
%             y(k)=aux(x(k));
%         end
%         plot(x,y)
%     end
%     b(i) = P(i-1)/P(i)*b(i-1);
%     b(i) = b(i)-br*l;
% end


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
%% -------------------Algoritmo de Chtcherbakov------------------------- %%
%normalización de potencia distribuída
P=P/Ptot;
clear br;
br=1.1;
ls = 0.795492119137066; % pérdida de acoplamiento obtenida de ajuste de curvas.
ls=0.8;
% ls=0.95
lm=1/1.01;
%-------- energía hacia la carga----%
% correccion = [0.00498654246545721,0.0173917101049116,0.0156910761023727,0.0134882951812974,0.0107544531644847,0.00763075701833112,0.00412958846942688,0.000757942896187622,-0.00202840856921822,-0.00473329344921942,-0.00644967660980973,-0.00788692650397274,-0.00827600057243928,-0.00825218757641784,-0.00733979259091601,-0.00666448259195723,-0.00578243355134984,-0.00575226995987445,-0.00559481514948685,-0.00669384250782859;]
% correccion2 =[-0.00150529963117793,-0.00497897790519619,0.00227767907607958,0.00804704709781173,0.0108225809233753,0.0117320502022904,0.0113138010645880,0.00888944018343983,0.00561760603509701,0.00205015051279949,-0.00125563685070079,-0.00416401629633608,-0.00588369543239822,-0.00689109206858135,-0.00663731165420330,-0.00631831913713530,-0.00564636479541522,-0.00570441485569356,-0.00562343805344107,-0.00673153212602525;]
% P=P+correccion;
P=0.995*P;
%----------------------------------%
cumP=cumsum(P);
%Vector de acoplamientos:
for i=[1:1:N]
    if (i==1)
        S2(i)=P(i); 
    else
        S2(i) = P(i) / (1 - cumP(i-1));
    end
end
% 

%  S2 =S2-correccion-correccion2;
b2=zeros(size(P)); %se inicializa vector de alturas  para el nuevo algoritmo

for i=[1:1:N-1]
    b2(i) = br* sqrt((ls^2*lm^(2*i-2)-ls*lm^(i-1)*S2(i))/(S2(i)*S2(i+1))) ;    
end
b2(N) = br* sqrt((ls^2*lm^(2*i-2)-ls*lm^(i-1)*S2(N))/(S2(N))) ;    

b0=ls*br/S2(1);
% % b2(1)=b2(1)*0.927861859705280;
% b2(1)=b2(1)*0.915;
% b2(2)=0.985915906645940*b2(2);% b2(2) = b2(2)-0.15;
% % b2(3)=b2(3)-0.3;
% 
% % b0=1.044558513037895*b2(1); %b0=b2(1)+0.5;
% b0=b2(1)*1.027732172638890;%b0=b2(1)+0.3;

  
  
%% Gráfico
figure(4);
step1 = -b+b0;
step2 = -b2+b0;
stairs([0:20]',[step1 step1(N);step2 step2(N)]');
legend('Método 1','Método 2');hold on;
plot(b0*ones(size(A)),'r');
title('Comparación de algoritmos');
ylabel('Altura de step (mm)');
xlabel('step');

%% --------------------------------------------------------------------- %%
%% --------------------------------------------------------------------- %%

% fout = fopen('completedist_results/steps_generic_cambiar.vbs','wt'); % open steps_generic_cambiar.vbs
fout = fopen('completedist_results/steps_cambiar.vbs','wt'); % open steps_generic_cambiar.vbs

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
fprintf(fout,'Set oDesign = oProject.SetActiveDesign("steps") \n');
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
fprintf(fout,'  Array("NAME:B", "Value:=", "%3.6gmm")))) ',b0);
fclose(fout);