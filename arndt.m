% Power distribution using Arndt algorithm
% Cristian Duguet S.
% 2009-12-21

clc
clear all;
close all;

%B: outputs
%A: inputs
%ST: Scattering matrix of the symmetrical T-junction
%SUT:  Scattering matrix of the unssymmetrical T-junction

%----------------------------Constants------------------------------------%
eps0=8.854187817e-12;
mu0=4*pi*10^-7;

%Permeability and Permitivity
epsr=1.0 ; mur=1; 
eps=eps0*epsr; 
mu=mu0*mur;
%----------------------Waveguide sizes------------------------------------%
% 3 zones of the unsymmetrical T-junction

%      br 
%  ___| |___| |___| |___| |___| |___| |___| |___| |___| |___| |___| |___| |___| |___| |___| |___| |___| |___| |___| |___| |___
%  b0    b1    b1    b3    b4    b5    b6    b7    b8    b9    b10   b11   b11   b13   b14   b15   b16   b17   b18   b19   b20|
%  _________________________________________-----------------------------------------------°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
%  
f0=10e9; %central frequency

%Waveguide width
a=15.799e-3;

%Height values

c=3e8/f0 /4;  %is NOT speed of light!!! is the distance between T's (lambda/4)
br=1.4e-3;
b0   = 11.0754e-3;
b(1) = 10.7765e-3;
b(2) = 10.3371e-3;
b(3) = 9.55055e-3;
b(4) = 8.79401e-3;
b(5) = 8.09139e-3;
b(6) = 7.42815e-3;
b(7) = 6.81862e-3;
b(8) = 6.25514e-3;
b(9) = 5.73560e-3;
b(10)= 5.24971e-3;
b(11)= 4.81309e-3;
b(12)= 4.45751e-3;
b(13)= 4.21883e-3;
b(14)= 4.08644e-3;
b(15)= 4.01240e-3;
b(16)= 3.83483e-3;
b(17)= 3.39190e-3;
b(18)= 2.60814e-3;
b(19)= 1.60244e-3;
b(20)= 0.64673e-3;   

% element counter 
w=0;
%frequency of analysis
for f=8.5e9:.5e9:11.5e9
    w=w+1;
    omega=2*pi*f;
    %wave number
    k= omega*sqrt(mu*eps);
    for N=1:20
        
        if N==1
            b1=b0;
        else
            b1=b(N-1);
        end
        b2=b(N);       
    %% -------------------ST MATRIX CALCULATION---------------------------%   
    ST(N,:,:)=zeros(3);
    
    for n = 0
        for p=0
            %------ We'll use TE_xn mode in I and TE_xp mode in III ------%
    
            %          |/\/\|
            %          |  p |
            % _________|    |__________
            %      /           \
            %    n \           / n
            % _____/___________\_______
            
            %Shear wave number, (for each mode)
            kx= pi/a;           %mode 1 ,x dependent, TE1_
            kq(1)= n*pi./b1;  %mode n , q dependent, TE1n
            kq(2)= kq(1);
            kq(3)= p*pi./br;  %mode p , q dependent, TE1n
            
            kc = sqrt(kx^2 + kq.^2);    %<-- 3-element vector (Shear wave number in I, II and III)
            
            
            %phase constant
            beta = sqrt(k^2 - kc.^2);
            beta(k^2 < kc.^2) = -1* beta(k^2 < kc.^2);   % 3-element vector
            
            
            %Impedance
            ZF = omega*mu./beta;                %3-element vector
            
            %Admitance
            Y = ZF.^(-1)                     %3-element vector
            
            
            %Diagonal Matrix e
            e=zeros(1,3);                %3-element vector
            e(1)=exp(1i*beta(1)*br);
            e(3)=exp(1i*beta(3)*b1);
            
            
            y1 = b1/2;
            z1 = br/2;
            
            %T Matrix, first element
            T(1) = 1/sqrt(y1*z1) * sin(beta(1)*2*z1)/sin(beta(3)*2*y1) * beta(1)/beta(3) * 1/sqrt(1+isequal(n,0)) *  1/sqrt(1+isequal(p,0)); %% me salta una duda, es beta(3)*2*z1?? o *y1???
            %multiplicado por la integral
            T(1)=T(1) * ((n*pi/2/y1)*sin(n*pi) - beta(3)*sin(beta(3)*2*y1)) / ((n*pi/2/y1)^2 - (beta(3))^2) ;
            %INTEGRAL OK
            
            %T Matrix, second element
            T(2) = T(1)*(-1)^p;  %% should it be the same power as T(4)?? p or n?
            
            %T Matrix, third element
            T(3) = 1/sqrt(y1*z1) * sin(beta(3)*2*y1)/sin(beta(1)*2*z1) * beta(3)/beta(1) * 1/sqrt(1+isequal(n,0)) *  1/sqrt(1+isequal(p,0));
            %multiplicado por la integral
            T(3)=T(3) * ((p*pi/2/z1)*sin(p*pi) - beta(1)*sin(beta(1)*2*z1)) / ((p*pi/2/z1)^2 - (beta(1))^2) ;
            
            %Matrix T, fourth element
            T(4) = T(3) * (-1)^p;
            
            
            %Simmetric T-junction Matrix
            STaux= [ -e(1)*Y(1)    Y(1)    -T(1)* Y(3);
                       Y(1)    -e(1)*Y(1)   T(2)* Y(3);
                     -T(3)*Y(1) T(4)*Y(1)   -e(3)*Y(3)];
            
            STaux = STaux \ [ (1/e(1))*Y(1)       -Y(1)          T(1)* Y(3);
                                  -Y(1)       (1/e(1))*Y(1)     -T(2)* Y(3);
                                T(3)*Y(1)      -T(4)*Y(1)     (1/e(3))*Y(3)];
            
            ST(N,:,:)=squeeze(ST(N,:,:))+STaux;
        end
    end
%     ST = ST/((n+1)*(p+1));
    STdb = mag2db(abs(ST));
    clear beta kq kc ZF Y
    
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    %% -------------------SH MATRIX CALCULATION---------------------------%

    SH(N,:,:)=zeros(2);
    for m =0
        for n =0
            %----- We'll use TE_xm mode in IIa and TE_xn mode in  V ------% 
            
            %    _____________________________
            %          \            / n
            %      m   /      ______\_________
            %    ______\_____|

            
            %Shear wave number, (for each mode)
            kx= pi/a;      %mode 1 ,x dependent, TE1_
            kq(1)= m*pi./b1;  %mode m , q dependent, TE1n       IN
            kq(2)= n*pi./b2;  %mode n , q dependent, TE1n      OUT

            kc = sqrt(kx^2 + kq.^2);    %<-- 2-element vector (Shear wave number in I and II)

            %phase constant
            beta = sqrt(k^2 - kc.^2);      
            beta(k^2 < kc.^2) = -1i* beta(k^2 < kc.^2);   % 2-element vector

      delay(N)=exp(1i*beta(2)*(c-br));       %Phase difference between 2 adjacent T-junctions in an array
            
            %Impedance
            ZF = omega*mu./beta;                %-element vector

            %Admitance
            Y = ZF.^(-0.5);                     %3-element vector

            %depends on b(i+1)
            K1 = 1/sqrt(1+isequal(m,0)) *  1/sqrt(1+isequal(n,0)) * 2/sqrt((2*y1)*(b2));

            %times the itegral
            if m==0 && n==0
              K1 = K1 * b2;
            elseif m==0
              K1 = K1 * sin(n*pi) * b2 / (n*pi);
            elseif n==0
              K1 = K1 * (sin(m*pi) - sin(m*pi*(b1-b2)/b1))* b1/(m*pi);
            else
              K1 = K1 * ( (m*pi/b1)*sin(m*pi)*cos(n*pi) - (n*pi/b2)*sin(n*pi)*cos(m*pi) - (m*pi/b1)*sin(m*pi*(b1-b2)/b1) ) / ( (m*pi/b1)^2 - (n*pi/b2)^2 );
            end;

            K2=transpose(K1); %transposed

            %diagonal matrices sqrtY 
            %Phase constant is the same between IIa and  I, and between  V and II:
            sqrtYm = Y(1);
            sqrtYn = Y(2);

            betam=beta(1);
            betan=beta(2);
            %-----------------------------%

            SHaux(1,1) =          ( sqrtYm * betam^(-1) * (K1*sqrtYm)^(-1) * sqrtYm +  K2 * betan^(-1) *sqrtYn )^(-1);
            SHaux(1,1) = SHaux(1,1) *(-sqrtYm * betam^(-1) * (K1*sqrtYm)^(-1) * sqrtYm +  K2 * betan^(-1) *sqrtYn );

            SHaux(1,2) =        2*( sqrtYm * betam^(-1) * (K1*sqrtYm)^(-1) * sqrtYm +  K2 * betan^(-1) *sqrtYn )^(-1) * sqrtYm * betam^(-1);

            SHaux(2,1) = 2*(  K1*sqrtYm  +  sqrtYn *(K2 * betan^(-1) *sqrtYn )^(-1) * sqrtYm *betam^(-1)  )^(-1) * sqrtYn;

            SHaux(2,2) =          (  K1*sqrtYm  +  sqrtYn *(K2 * betan^(-1) *sqrtYn )^(-1) * sqrtYm *betam^(-1)  )^(-1);
            SHaux(2,2) = SHaux(2,2) *( -K1*sqrtYm  +  sqrtYn *(K2 * betan^(-1) *sqrtYn )^(-1) * sqrtYm *betam^(-1)  );

%             if m==0 && n==0
            SH(N,:,:)= squeeze(SH(N,:,:))+SHaux;
%             else 
%                 SH=SH + SHaux/10;
%             end
            %------------------------------%
        end
    end
%     SH=SH/((m+1)*(n+1));
    SHdb= mag2db(abs(SH));
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    
    %% -------------TOTAL MODAL SCATTERING MATRIX-------------------------%


    M1 = (ones(1) - ST(N,2,2)* SH(N,1,1))^(-1);
    M2 = SH(N,1,1)*M1;
    M3 = SH(N,2,1)*M1;
    M4 = ST(N,2,2)*SH(N,1,2);
    %-------------------------------%
    SUT(N,1,1) = ST(N,1,1) + ST(N,1,2)*M2*ST(N,2,1);
    SUT(N,1,2) = ST(N,1,2)*SH(N,1,2) + ST(N,1,2)*M2*M4;
    SUT(N,1,3) = ST(N,1,3) + ST(N,1,2)*M2*ST(N,2,3);
    SUT(N,2,1) = M3*ST(N,2,1);
    SUT(N,2,2) = SH(N,2,2) +M3*M4;
    SUT(N,2,3) = M3*ST(N,2,3);
    SUT(N,3,1) = ST(N,3,1) + ST(N,3,2)*M2*ST(N,2,1);
    SUT(N,3,2) = ST(N,3,2)*SH(N,1,2) + ST(N,3,2)*M2*M4;
    SUT(N,3,3) = ST(N,3,3) + ST(N,3,2)*M2*ST(N,2,3);
    end %end of T elements
    
    %S-parameters of the array
    S(1,w) = SUT(1,3,1);
    aux=cumprod(SUT(:,2,1).*delay(:)); %cumulative product of intermediate T's and delays
    for N=2:20
        S(N,w) = aux(N-1) * SUT(N,3,1);
    end

end %end of frequency sweep

%-------------------------------------------------------------------------%
%Distribución de Taylor Teórica,  SLL=30dB, n_=7
A=[40.0173   39.6844   39.0521   38.1584   37.0171   35.6173   33.9545   32.0634   30.0176   27.8902   25.7070   23.4396   21.0569   18.6039   16.2431   14.2140   12.7217   11.8276   11.4199   11.2941];
% Dist de Potencias de Taylor
P=A.^2;
P=P/sum(P);
%-------------------------------------------------------------------------%
Scuad=abs(S).^2;


freq = 8.5e9:.5e9:11.5e9;

plot(1:20,Scuad(:,4),1:20,P); legend('Arndt','Goal');

figure(2)
plot(1:20,rad2deg(angle(S(:,4))))





% s11 = [-29.586039 -29.317956 -32.413265 -37.012861 -37.209029 -31.812850 -27.553865 -24.317514 -21.684441 -19.432421 -17.434651 -15.613565 -13.919435 -12.320532 -10.798725 -9.347511 -7.970104];
%     
% % figure(1)
% plot(freq,mag2db(abs(SUT(:,1,1))),freq,s11);legend('S_1_1 Arndt','S_1_1 simulado'); xlabel('freq [GHz]'); ylabel('[dB]'); title('Reflexión')
% 
% s21 = [-2.821004 -2.905269 -2.876396 -2.813972 -2.751084 -2.689417 -2.626774 -2.562317 -2.496345 -2.429624 -2.363711 -2.300775 -2.243744 -2.196838 -2.166497 -2.162548 -2.199255];
% s21 = [-2.913054 -2.898496 -2.857597 -2.806131 -2.749507 -2.689450 -2.626774 -2.562215 -2.496893 -2.432625 -2.372135 -2.319198 -2.278679];
% s21 = [-3.020781 -2.912594 -2.858067 -2.805238 -2.748866 -2.689137 -2.626763 -2.562310 -2.496294 -2.429427 -2.362909 -2.298714 -2.239919 -2.191087 -2.158802 -2.152509 -2.185866];
% s21 = [-2.963696 -2.913711 -2.861366 -2.806502 -2.749066 -2.689106 -2.626774 -2.562349 -2.496294 -2.429354 -2.362731 -2.298332 -2.239154 -2.189810 -2.157253 -2.151721 -2.187890]; 
% 
% s31 = [-3.274481 -3.145905 -3.175725 -3.228926 -3.290462 -3.363315 -3.447845 -3.546469 -3.663499 -3.802332 -3.968078 -4.168246 -4.413719 -4.720077 -5.109236 -5.611093 -6.264045]; 
% s31 = [-3.073984 -3.118983 -3.169236 -3.225822 -3.290011 -3.363361 -3.447845 -3.546030 -3.661310 -3.798214 -3.962827 -4.163347 -4.410854 -4.720376 -5.112358 -5.614577 -6.264282];
% 
% figure(2)
% plot(freq,mag2db(abs(SUT(:,2,1))),freq,mag2db(abs(SUT(:,3,1))),freq,s21,'.-',freq,s31,'.-');
% legend('S_2_1 Arndt','S_3_1 Arndt', 'S_2_1 simulado','S_3_1 simulado'); 
% title('Parametros S obtenidos con Algoritmo Arndt'); xlabel('freq [GHz]'); ylabel('[dB]');

% figure(3)
% plot(freq,rad2deg(angle(SUT(:,2,1))),freq,rad2deg(angle(SUT(:,3,1))));legend('S_2_1','S_3_1');


