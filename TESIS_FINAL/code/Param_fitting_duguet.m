%Cristian Duguet

%Ajuste Paramï¿½trico del la curva multidimensional de acoplamiento de un
%divisor de potencia.

%Modificaciones:
%2009-oct-12: C.D. - Creaciï¿½n
%2009-oct-23: C.D. - ObtenciÃ³n de parÃ¡metro de reflexiÃ³n
%-------------------------------------------------------------------------%
clc
clear all
close all
%% ----------------------- IMPORTACION DE DATOS --------------------------%
%Datos de simulación importados  en step_parametric/Sparam_generic.csv (para div cuadrado)
newData1 = importdata('step_parametric/step_generic_S_br1.4.csv');           
% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

%Datos de simulaciï¿½n importados  en S_generic.mat (para div cuadrado)
% open('S_generic.mat');
%Datos de simulaciï¿½n importados  en S_sloping.mat (para div inclinado)
%open('S_sloping.mat');
% data=ans.data;
%-------------------------------------------------------------------------%
clear newData1 vars textdata

%Creaciónn de la matrix 4-dimensional: las columnas dependientes de la 
%frecuencia en la matriz importada estan ordenados en br, dentro en b1, 
%dentro en b.
for ibr=1; %br=1.4
for ib1=1:12; %b1=1,...,12
for ib=ib1:12; %b=2,...,13
    S_generic(ib1,ib,ibr,:)=data(:,157+((ibr-1)*45)+((27*ib1-ib1^2)/2-12)+(ib-ib1));
end
end
end

%Creaciónn de la matriz de inserciones. Acoplamiento hacia el siguiente
%divisor, out1.
for ibr=1; %br=1.4
for ib1=1:12; %b1=1,...,12
for ib=ib1:12; %b=2,...,13
    S_ins(ib1,ib,ibr,:)=data(:,79+((ibr-1)*78)+((27*ib1-ib1^2)/2-12)+(ib-ib1));
end
end
end

%Matriz de reflexiones
for ibr=1; %br=1.4
for ib1=1:12; %b1=1,...,12
for ib=ib1:12; %b=2,...,13
    S_ref(ib1,ib,ibr,:)=data(:,1+((ibr-1)*78)+((27*ib1-ib1^2)/2-12)+(ib-ib1));
end
end
end


%% Independencia con respecto de b
%Para ver la independencia con respecto de b, y si es bueno considerar que
%b=b1+br. Para esto se utilizarï¿½ el criterio de mï¿½ximo  acoplamiento (S)


                    %------------ Tabla ------------%
for ibr = 1
    fprintf('Valor de S máximos, para cada b_1, b_r = %d mm , f =10GHz \n',ibr*1.4)    
    [maxi ind] = max(S_generic(:,:,ibr,16),[],2); %los máximos para cada b1 fixed
    for ib1=1:9
        mini(ib1) = min(  S_generic(ib1,S_generic(ib1,:,ibr,16) > 0,ibr,16)); %los mínimos distintos de 0
        valor_restr(ib1) = S_generic(ib1,min(ibr+ib1+1,9),ibr,16);
    end
    
    fprintf('b1 \t & \t b (S máx) \t & \t Variación máxima(\\%%) \t & \t Variación en b=b1+br (\\%%) \\\\ \\hline \n');
    for ib1=1:9
        fprintf('%g \t & \t %g \t & \t %3.3f \t & \t %3.4f  \\\\ \\hline \n',ib1, ind(ib1), (maxi(ib1)-mini(ib1))/maxi(ib1), (maxi(ib1)-valor_restr(ib1))/maxi(ib1))
    end
end

clear mini maxi

%Entonces, de lo anterior se ve que la respuesta de S tiende a un máximo
%cuando b tiende a la cercanía de b1+br, sin embargo esta respuesta sigue
%siendo muy plana. Se recomienda, para simplificar el modelo a ajustar, que
%se considere la restricción b = b1+br o cercana, y se consideren sólo las
%variables b1,br,f en el ajuste paramétrico de curvas, de acuerdo al modelo
%experimental.
%Esto ademas nos facilita el ajujste de curvas del modelo, ya que
%considerando b, existia la condición b1<b, entonces había que realizar un
%ajuste en una matriz 4-dimensional no "cuadrada". 

%-------------------------------------------------------------------------%
%Se modela la matriz de datos restringida como:

for ibr =1 %for br = 1.4mm
    for ib1 = 1:12 %for b1 =1mm,...,12mm
        S_generic_restr(ib1,ibr,:) = S_generic(ib1,min(ib1+ibr+1,9),ibr,:);
    end
end

%y la matriz de reflexiones:
for ibr =1 %for br = 1.4mm
    for ib1 = 1:12 %for b1 =1mm,...,12mm
        S_ins_restr(ib1,ibr,:) = S_ins(ib1,min(ib1+ibr+1,9),ibr,:);
    end
end

% clear S_generic;
% clear S_reflex;

%-------------------------------------------------------------------------%

%% --------------------------------------------------------------------- %%
% --------------------OBTENCION DE PARAMETROS PARA EL-------------------- %
% ------------------ALGORITMO DUGUET (l y alpha)------------------------- %
% ----------------------------------------------------------------------- %

%Recordar que las variables son dependientes de : (b1,b,br,f)

%Potencia de salida a radiador:
P_r = S_generic.^2;
%Potencia hacia la siguiente guÃ­a:
P_b1 = S_ins.^2;
%Potencia reflejada
P_reflex = S_ref.^2;

P_r+P_b1+P_reflex

%factor alpha de atenuaciÃ³n y reflexiones:
alpha = 1./(P_r + P_b1);

datos = squeeze(alpha(:,:,1,16));
[a,b] = ndgrid(1:12,2:13);
% [a,b] = ndgrid(3:7, 8.5:.1:11.5);
figure(1)
surface(a,b,datos);xlabel('b1');ylabel('b');hold on;
% datos = squeeze(alpha(:,:,1,1));
% surface(a,b,datos);xlabel('b1');ylabel('b');hold on;
% datos = squeeze(alpha(:,:,1,31));
% surface(a,b,datos);xlabel('b1');ylabel('b');hold on;


%se observa que alpha depende de la forma de b-b1, entonces se restringe el
%parámetro a b-b1, independiente de la frecuencia, y con máximo en b-b1=br


for i=0:5;
alpha_restr(1+2*i,:,:) = alpha(7-i,7+i,:,:);
alpha_restr(2+2*i,:,:) = alpha(6-i,7+i,:,:);
end
% alpha_restr(1,:,:) = alpha(7,7,:,:);
% alpha_restr(2,:,:) = alpha(6,7,:,:);
% alpha_restr(3,:,:) = alpha(6,8,:,:);
% alpha_restr(4,:,:) = alpha(5,8,:,:);
% alpha_restr(5,:,:) = alpha(5,9,:,:);
% alpha_restr(6,:,:) = alpha(4,9,:,:);
% alpha_restr(7,:,:) = alpha(4,10,:,:);
% alpha_restr(8,:,:) = alpha(3,10,:,:);
% alpha_restr(9,:,:) = alpha(3,11,:,:);
% alpha_restr(10,:,:) = alpha(2,11,:,:);
% alpha_restr(11,:,:) = alpha(2,12,:,:);
% alpha_restr(12,:,:) = alpha(1,12,:,:);


alpha_coef = polyfit ([1:12]',alpha_restr(:,16),3);
x=1:12;
figure(2)
plot(x,alpha_coef(1)*x.^3+alpha_coef(2)*x.^2+alpha_coef(3)*x+alpha_coef(4),x,alpha_restr(:,16))

figure(3)
x=b-a;
surface(a,b,datos);xlabel('b1');ylabel('b');hold on;
surface(a,b,alpha_coef(1)*x.^3+alpha_coef(2)*x.^2+alpha_coef(3)*x+alpha_coef(4));xlabel('b1');ylabel('b');

%-------------------------------------------------------------------------%
%%Parámetro l
[ib1,ib,ibr,i_f] = ndgrid(1:12,2:13,1,8.5:.1:11.5);
br=1.4;
l = (P_r./P_b1).*ib1/br;
b=isnan(l);
l(b)=0;
% ---------l vs frecuencia, br----------%
figure(4)
datos = squeeze(l(3,:,1,:)); [a,b] = ndgrid(2:13,8.5:.1:11.5);
surface(a,b,datos);xlabel('b');ylabel('f');hold on;

% Claramente, l depende linealmente de br, sin embargo también varía con la
% frecuencia. Para b_r mas pequeños (\approx 1.4mm) , l es más plano en
% frecuencia.

%---------l vs b1,b          ----------%
figure(5)
datos = squeeze(l(:,:,1,16));
[a,b] = ndgrid(1:12,2:13);
surface(a,b,datos);xlabel('b1');ylabel('b');hold on;

% Ajuste de parámetros
l_coef = polyfit ([1:12]',l(:,12,1,16),1);



