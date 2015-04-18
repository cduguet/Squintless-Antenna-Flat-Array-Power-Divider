%Cristian Duguet

%Ajuste Paramétrico del la curva multidimensional de acoplamiento de un
%divisor de potencia.

%Modificaciones:
%2009-oct-12: C.D. - Creación
%2009-oct-23: C.D. - ObtenciÃ³n de parÃ¡metro de reflexiÃ³n
%-------------------------------------------------------------------------%
clc
clear all
close all
%% ----------------------- IMPORTACION DE DATOS --------------------------%
%Datos de simulación importados  en step_parametric/Sparam_generic.csv (para div cuadrado)
newData1 = importdata('step_parametric/step_slope_091026.csv');           
% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

%Datos de simulación importados  en S_generic.mat (para div cuadrado)
% open('S_generic.mat');
%Datos de simulación importados  en S_sloping.mat (para div inclinado)
%open('S_sloping.mat');
% data=ans.data;
%-------------------------------------------------------------------------%
clear newData1;
clear vars;

%Creación de la matrix 4-dimensional: las columnas dependientes de la 
%frecuencia en la matriz importada estan ordenados en br, dentro en b1, 
%dentro en b.
for ibr=1:5; %br=3,4,5,6,7
for ib1=1:9; %b1=1,...,9
for ib=ib1:9; %b=2,...,10
    S_generic(ib1,ib,ibr,:)=data(:,226+((ibr-1)*45)+((21*ib1-ib1^2)/2-9)+(ib-ib1));
end
end
end

%Creación de la matriz de inserciones. Son exactamente los mismos
%parómetros de acoplamiento, pero asociados a la salida "out1" que es la
%salida que conecta al siguiente divisor. 
for ibr=1:5; %br=3,4,5,6,7
for ib1=1:9; %b1=1,...,9
for ib=ib1:9; %b=2,...,10
    S_ins(ib1,ib,ibr,:)=data(:,1+((ibr-1)*45)+((21*ib1-ib1^2)/2-9)+(ib-ib1));
end
end
end
%-------------------------------------------------------------------------%
%Para ver la independencia con respecto de b, y si es bueno considerar que
%b=b1+br. Para esto se utilizaró el criterio de móximo  acoplamiento (S)

for ibr = 1:5
    fprintf('Valor de S móximos, para cada b_1, b_r = %d mm , f =10GHz \n',ibr+2)    
    [maxi ind] = max(S_generic(:,:,ibr,16),[],2); %los móximos para cada b1 fixed
    for ib1=1:9
        mini(ib1) = min(  S_generic(ib1,S_generic(ib1,:,ibr,16) > 0,ibr,16)); %los mónimos distintos de 0
        valor_restr(ib1) = S_generic(ib1,min(ibr+ib1+1,9),ibr,16);
    end
    
    fprintf('b1 \t & \t b (S móx) \t & \t Variación móxima(\\%%) \t & \t Variación en b=b1+br (\\%%) \\\\ \\hline \n');
    for ib1=1:9
        fprintf('%g \t & \t %g \t & \t %3.3f \t & \t %3.4f  \\\\ \\hline \n',ib1, ind(ib1), (maxi(ib1)-mini(ib1))/maxi(ib1), (maxi(ib1)-valor_restr(ib1))/maxi(ib1))
    end
end

figure(5)
[a,b] = ndgrid(2:10,1:9);
data=S_generic(:,:,3,16);
plot(a,data'); title('Variación del acoplamiento, con b');
xlabel('b');
for ib1=1:9 
    legend(sprintf('b1=%g mm',ib1));
end

%Entonces, de lo anterior se ve que la respuesta de S tiende a un móximo
%cuando b tiende a la cercanóa de b1+br, sin embargo esta respuesta sigue
%siendo muy plana. Se recomienda, para simplificar el modelo a ajustar, que
%se considere la restricción b = b1+br o cercana, y se consideren sólo las
%variables b1,br,f en el ajuste paramótrico de curvas, de acuerdo al modelo
%experimental.
%Esto ademas nos facilita el ajujste de curvas del modelo, ya que
%considerando b, existia la condición b1<b, entonces habóa que realizar un
%ajuste en una matriz 4-dimensional no "cuadrada". 

%-------------------------------------------------------------------------%
%Se modela la matriz de datos restringida como:

for ibr =1:5 %for br = 3mm...5mm
    for ib1 = 1:9 %for b1 =1mm,...,9mm
        S_generic_restr(ib1,ibr,:) = S_generic(ib1,min(ib1+ibr+1,9),ibr,:);
    end
end

%y la matriz de reflexiones:
for ibr =1:5 %for br = 3mm...5mm
    for ib1 = 1:9 %for b1 =1mm,...,9mm
        S_ins_restr(ib1,ibr,:) = S_ins(ib1,min(ib1+ibr+1,9),ibr,:);
    end
end

% clear S_generic;
% clear S_reflex;

%-------------------------------------------------------------------------%
%% Ajuste de parómetros del modelo a los datos
%Para esto se utiliza el mótodo de ajuste por mónimos cuadrados con
%lsqcurvefitting. Las variables independientes son 3: b1, br, f
%Los parómetros son 2.

% La ecuacion empórica tiene la forma matemótica de: 
% S^2 = br/(b) * alpha * (10/f^beta)
%   con b=b1+br:
% S^2 = br/(b1+br) * alpha * (10/f^beta)

%Se crean las matrices de variables independientes del modelo, usando
%ndgrid
[xdata.b1,xdata.br,xdata.f] =   ndgrid(1:1:9, 3:1:7, 8.5:.1:11.5);

%Los datos obtenidos de la simulación. Datos a los que se ajusta la curva.
ydata=S_generic_restr.^2;

%Valores iniciales de parómetros para la iteración. Pero como se trata de
%la modelacion de una estructura independiente de la frecuencia, debemos
%encontrar el ajuste de minimos cuadrados para el modelo, sin la variable
%f. Entonces el parametro potencia de f, desaparece.
x0(1) = 0.3;
%x0(2) = .5;

%DEFINICION DEL MODELO, función de parómetros (x) y puntos de variables
%independientes (x_var). Cabe destacar que los nombres son derivados de la
%nomenclatura usada en MATLAB.
Fun = @(x,x_var)x_var.br./(x_var.br+x_var.b1) *x(1)%.*10./x_var.f.^x(2);


%AJUSTE CON MINIMOS CUADRADOS
[x, resnorm] = lsqcurvefit(Fun, x0, xdata, ydata)

%-------------------------------------------------------------------------%
%GRóFICO 1 : modelo y datos, vs  b1 , f ---- br=3mm
figure(1)
[a,b]=ndgrid(1:9, 8.5:.1:11.5);
% x= [1/1.4 0]
datos = squeeze(ydata(:,1,:));
modelo = Fun(x,xdata);
modelo=squeeze(modelo(:,1,:));
surface(a,b,modelo);hold on;surface(a,b,datos); hold off;
% plot3(a,b,modelo,a,b,datos);
xlabel('b1');ylabel('f');

%-------------------------------------------------------------------------%
%GRóFICO 1 : modelo y datos, vs  b1 , br ------ f=10GHz
figure(2)
[a,b]=ndgrid(1:9, 3:7);
datos_br = squeeze(ydata(:,:,16));
modelo = Fun(x,xdata);
modelo=squeeze(modelo(:,:,16));
surface(a,b,modelo);hold on;surface(a,b,datos_br);xlabel('b1');ylabel('br'); hold off;

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%% Ajuste Paramótrico de Inserción:
%Se trataró de ajustar paramótricamente curvas que representen a la
%Potencia transmitida a travós del divisor hacia el divisor siguiente, es 
%decir, considerando las reflexiones, o la pórdida de inserción.
%Para esto se utiliza el mótodo de ajuste por mónimos cuadrados con
%lsqcurvefitting. Las variables independientes son 3: b1, br, f
%Los parómetros son 2.

% La ecuacion empórica tiene la forma matemótica de: 
% S^2 = b1/(b) * gamma * (10/f^kappa)
%   con b=b1+br:
% S^2 = b1/(b1+br) * gamma * (10/f^kappa)

%Los datos obtenidos de la simulación. Datos a los que se ajusta la curva.
ydata_ins=S_ins_restr.^2;

%Valores iniciales de parómetros para la iteración. Pero como se trata de
%la modelacion de una estructura independiente de la frecuencia, debemos
%encontrar el ajuste de minimos cuadrados para el modelo, sin la variable
%f. Entonces el parametro potencia de f, desaparece.
x0_ins(1) = 0.3;
%x0_ins(2) = .5;

%Definición del modelo de inserción de potencia al divisor siguiente, 
%función de parómetros (x) y puntos de variables
%independientes (x_var). Cabe destacar que los nombres son derivados de la
%nomenclatura usada en MATLAB.
Fun_ins = @(x,x_var)x_var.b1./(x_var.br+x_var.b1) *x(1)%.*10./x_var.f.^x(2);


%AJUSTE CON MINIMOS CUADRADOS
[x_ins, resnorm_ins] = lsqcurvefit(Fun_ins, x0_ins, xdata, ydata_ins)

%-------------------------------------------------------------------------%
%GRóFICO 1 : modelo y datos, vs  b1 , f ---- br=3mm
figure(3)
[a,b]=ndgrid(1:9, 8.5:.1:11.5);
datos = squeeze(ydata_ins(:,1,:));
modelo_ins = Fun_ins(x_ins,xdata);
modelo_ins=squeeze(modelo_ins(:,1,:));
plot3(a,b,modelo_ins,a,b,datos);xlabel('b1');ylabel('f');

%-------------------------------------------------------------------------%
%GRóFICO 1 : modelo y datos, vs  b1 , br ------ f=10GHz
figure(4)
[a,b]=ndgrid(1:9, 3:7);
datos_br = squeeze(ydata_ins(:,:,16));
modelo_ins = Fun_ins(x_ins,xdata);
modelo_ins=squeeze(modelo_ins(:,:,16));
surface(a,b,modelo_ins);hold on;surface(a,b,datos_br);xlabel('b1');ylabel('br'); hold off;

