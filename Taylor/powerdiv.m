clear all;
close all;
clear figures;
clc;
% Cristian Duguet Sáez. 2009-09-12
% Algoritmo para el cálclo de altura de steps en un divisor de potencia squintless
% 20 elementos
% Distribución, coseno cuadrado sobre pedestal
% Recordar que se realiza la distribución de potencia para la mitad del
% arreglo (arrelgo final tiene 40 elementos)
el=20;
%freq: 10GHz;
%lambda = 30 %mm;

%---dist coseno cuadrado--------%
%  x=0:(pi/2)/(el-1):pi/2;
% A=0.1+0.9.*(cos(x)).^2;
%  A=-x+pi/2;
%-------------------------------%

%%%%%%%%%%%%-----VALORES IMPORTADOS DE DIST. TAYLOR-----------%%%%%%%%%%%%%
%n=20, SLL=30dB
A=[38.8482   38.5671   38.0105   37.1863   36.1113   34.7999   33.2792   31.5677   29.7022  27.7007   25.6106   23.4418   21.2570   19.0430   16.8973   14.7367   12.7930   10.6124    10.4550   19.3801] ;

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

% for i=[1:1:N-1]
%     a(i)=P(i+1)/P(i);
%     c(i)= (Ptot - sum(P(1:i+1)))/(Ptot - sum(P(1:i)));
%     if a(i) < c(i)
%         sprintf('ERROR: La distribución de potencias decae muy rápido, o tiene variaciones muy rápidas')
%         return
%     end
% end

b0=13; %mm 
br=1.4; %mm

b1min=br * (Ptot-P(1))/P(1);
b1max=br * (Ptot-P(1))/(P(1)-P(N));
b(1)=b1min;

% VERSION ANTIGUA
% % b(1)= 12.5; %mm (altura de la guía en el primer escalón)
% % 
% % %--------------------------límites de br----------------------------------%
% % brmin=P(1)*b(1)/(Ptot - P(1)) - P(N)*b(1)/(Ptot-P(1))
% % brmax=P(1)*b(1)/(Ptot - P(1))
% % br = brmax %mm (altura de la guía que alimenta a cada radiador)
% % %-------------------------------------------------------------------------%

%-------------------- DETERMINACION DE b() -------------------------------%
for i=2:N
    b(i) = P(i-1)/P(i)*b(i-1);
    b(i) = b(i)-br;
end
%-------------------------------------------------------------------------%

%---------------------- NO   CRECIMIENTO DE B ----------------------------%
figure(1)
limit(1)=Ptot*br/(br+b(1));
for i=[2:N]
    b(i) = min([b(i) b(i-1)]);
    limit(i)=  (Ptot- sum(P(1:i-1))) * br/(br+b(i));
end
    
%Cómo se ve afectada la distribución por esta restricción:
figure(1)
semilogy([1:20],P,'b',[1:20],limit,'r');
title('Distribucion de Potencia'); legend('Distrib. de Taylor','Sin decrecimiento de steps');
%-------------------------------------------------------------------------%
clear a,i;
%% -------------------Algoritmo de Chtcherbakov------------------------- %%
%normalización de potencia distribuída
P=P/Ptot;
clear br;
br=1.4;

cumP=cumsum(P);
%Vector de acoplamientos:
for i=[1:1:N]
    if (i==1)
        S2(i)=P(i); 
    else
        S2(i) = P(i) / (1 - cumP(i-1));
    end
end

b2=zeros(size(P)); %se inicializa vector de alturas  para el nuevo algoritmo

for i=[1:1:N-1]
    b2(i) = br* sqrt((1-S2(i))/(S2(i)*S2(i+1))) ;    
end
  b2(N) = br* sqrt((1-S2(N))/(S2(N))) ;    

  b0=br/S2(1)
  
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
% 
%   Columns 1 through 9
% 
%    12.5000   11.3229   10.2969    9.3984    8.6063    7.9072    7.2863    6.7378    6.2508
% 
%   Columns 10 through 18
% 
%     5.8267    5.4565    5.1529    4.9066    4.7538    4.6778    4.6778    4.6778    4.6778
% 
%   Columns 19 through 20
% 
%     4.6778    0.0134
