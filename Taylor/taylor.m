%taylor.m
%distribucion de taylor para 6 elementos
clear all
 close all
clc

n=20;
SLL=30;%dB
b=10^(SLL/20);
A= acosh(31.6228) / pi;

for i=[1:n-1]
    UN(i) = n * sqrt(A^2 + (i - 0.5)^2) / sqrt(A^2 + (n - 0.5)^2);
end
clear i;
%Ceros en el espacio theta
% zeros= rad2deg(asin(UN/40/0.7))
B0=1;
%B coefficients
for i=[1:n-1]
    aux=1- i^2* [1:n-1].^(-2);%para calcular los coeficientes B
    aux(i)=1; % aux(i) es elemento absorbente y no debe estar en prod
    
    B(i) = 2* ( (-1)^i * prod(1- i^2* UN.^(-2)) ) / ( -2 * prod(aux));
end
clear i;

%-------------------------------------------------------------------------%
% No requiere normalizacion para el u-space pattern, pero si para calcular
% la distribucion de apertura
% B0=B0 /(B0 + sum(B));
% for k=[1:n-1]
%     B(k) = B(k) *B0;
% end
%-------------------------------------------------------------------------%

%-------------------------U-space pattern---------------------------------%
U=[0:0.001:10];
Taylor= B0*sin(pi*U)./(pi*U);
for i=1:n-1
    Taylor=Taylor + 0.5* B(i) *( sin(pi.*(U-i))./(pi.*(U-i)) + sin(pi.*(U+i))./(pi.*(U+i))); 
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
TaylordB=mag2db(abs(Taylor));
[tayHPBW,i] = min(abs(TaylordB+3));
tayHPBW=U(i); %Valor de HPBW variable U

%Factor de ancho de haz, considerando qque el ancho de haz nominal (dist
%uniforme) es  asin(0.4429lambda/a)

HPBWfactor= tayHPBW/0.4429
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(U,mag2db(abs(Taylor))); title('Patrón de Radiación en el Espacio-U');%hold on
xlabel('U (k sin \theta)');
ax=axis;
ax(3:4) = [-60 10];
axis(ax);
grid on;


%% ------------------ Aperture Distribution----------------------------- %%


x=0:0.001:0.5;
E = B0;
for i=[1:n-1]
    E = E + B(i).* cos(2*pi*i.*x);
end
figure(2)
plot(x,E); %hold on

%Entonces, se realiza un rough estimate de la potencia en cada elemento del 
%arreglo lineal, integrando la potencia de la disttribución continua en 20 elementos.          

%A: Voltaje en cada elemento

div=floor(size(E,2)/n); %ojo que toma un numero de puntos multiplo de n-1
for i=[1:n]
    A(i)= sum(E((i-1)*div+1:i*div+1));
end;
stem(A/sum(A)); title('Distribución normalizada de elementos del Arreglo, n=20, SLL=30dB')
xlabel('Elemento'); grid on;
    
%-------------------------------------------------------------------------%