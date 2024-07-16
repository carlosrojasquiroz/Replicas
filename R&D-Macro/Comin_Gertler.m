var Y Yc Yk Pc Lc Lk L K Psi Kc Kk U Mc Mk D Pk Kc Kk C R Muc Muk Zc Zk 
    Lambdac Lambdak Ac Ak Sc Sk Hc Hk Vc Vk PMc PMk PI PK Y Nc Nk Wr
parameters beta xi
beta=0.96;
xi=1;
delta=0.025*4;
eta=1;
    
%Consumidores

%Los consumidores consumen, trabajan e invierten en la producción y
%adopción de I+D.

%Condición intertemporal del consumo
C^(-1)=beta*C(+1)^(-1)*(D+PK(+1))/PK; %C
%Condición intratemporal del consumo
Wr=Muw*L^xi*C; %Wr
%Condición de arbitraje entre capital y préstamos
R(+1)=(D+PK(+1))/PK; %R
%Precios de bienes al consumidor
Pc=1; %Pc
%Ecuación de evolución del capital
K(+1)=(1-delta*U^eta)*K+Yk; %K

%Productores

%PMg del trabajo del sector de bienes al consumidor
(1-alpha)*(1-gamma)*Pc*Yc/Lc=Muc*Muw*L^xi*C; %Lc
%PMg del trabajo del sector de bienes de capital
(1-alpha)*(1-gamma)*Pk*Yk/Lk=Muk*Muw*L^xi*C; %Lk
%Oferta laboral total
L=Lc+Lk; %L
%Definición de capital del sector de bienes al consumidor
Kc=K/L*Lc; %Kc
%Definición de capital del sector de bienes de capital
Kk=K/L*Lk; %Kk
%Productividad marginal del capital del sector de bienes de capital
alpha*(1-gamma)*Pc*Yc/Kc=Muc*(D+delta*U^eta*PK); %D
%Productividad marginal del capital del sector de bienes al consumidor
alpha*(1-gamma)*Pk*Yk/Kk=Muk*(D+delta*U^eta*PK); %Pk
%Ratio de uso de capital
alpha*(1-gamma)*Pc*Yc/U =Muc*delta*eta*U^(eta-1)*PK*Kc; %U
%Precio relativo del capital
alpha*(1-gamma)*Pk*Yk/U =Muk*delta*eta*U^(eta-1)*PK*Kk; %PK
%Locación de bienes intermedios al sector de bienes al consumidor
gamma*Pc*Yc/Mc=Muc*PMc; %Mc
%Locación de bienes intermedios al sector de bienes de capital
gamma*Pk*Yk/Mk=Muk*PMk; %Mk
%Factor de escala
Psi=PI*K; %Psi
%Margen de ganancias del sector de bienes al consumidor
Muc=1/(Nc^mu); %Muc, inicialmente mu=1
%Margen de ganancias del sector de bienes de capital
Muk=1/(Nk^mu); %Muk, inicialmente mu=1
%Número de firmas de bienes finales al consumidor (entrada y salida)
(Muc-1)/Muc*Nc^(-Muc)*Yc=bx*Psi; %Nc
%Número de firmas de bienes finales de capital (entrada y salida)
(Muk-1)/Muk*Nk^(-Muk)*Yk=bx*Psi; %Nk
%Stock de innovación en el sector de bienes al consumidor
Zc(+1)/Zc=Chic*(Sc/Psi)^rho+phi; %Zc
%Stock de innovación en el sector de bienes de capital
Zk(+1)/Zk=Chik*(Sk/Psi)^rho+phi; %Zk


Y=C+Pk*Yk+G+(Sk+(Zk-Ak)*Hk)+(Sc+(Zc-Ac)*Hc); %Y
Yc=Nc^(Muc-1)*((U*K/L)^alpha*Lc)^(1-gamma)*Mc^gamma; %Yc
Yk=Nk^(Muk-1)*((U*K/L)^alpha*Lk)^(1-gamma)*Mk^gamma; %Yk




Lambdac=1/lambda*(Ac/Psi)*Hc^lambda; %Lambdac, con lambda<1
Lambdak=1/lambda*(Ak/Psi)*Hk^lambda; %Lambdak, con lambda<1  
Ac(+1)/Ac=Lambdac*phi*(Zc/Ac-1)+phi; %Ac
Ak(+1)/Ak=Lambdak*phi*(Zk/Ak-1)+phi; %Ak
R(+1)^(-1)*phi*Jc(+1)*(Zc(+1)-phi*Zc)=Sc; %Sc
R(+1)^(-1)*phi*Jk(+1)*(Zk(+1)-phi*Zk)=Sk; %Sk
Ac/Psi*(Ac/Psi)*Hc^(lambda-1)*R(+1)^(-1)*phi*(Vc(+1)-Jc(+1))=1; %Hc  
Ak/Psi*(Ak/Psi)*Hk^(lambda-1)*R(+1)^(-1)*phi*(Vk(+1)-Jk(+1))=1; %Hk
Vc=(Muc^(-1)*(1-1/nu)*gamma*Pc*Yc/Ac)+R(+1)^(-1)*phi*Vc(+1); %Vc
Vk=(Muk^(-1)*(1-1/nu)*gamma*Pk*Yk/Ak)+R(+1)^(-1)*phi*Vk(+1); %Vk
Jc=-Hc+R(+1)^(-1)*phi*(Jc(+1)+Lambdac*(Vc(+1)-Jc(+1)); %Jc
Jk=-Hk+R(+1)^(-1)*phi*(Jk(+1)+Lambdak*(Vk(+1)-Jk(+1)); %Jk
PMc=nu*(Ac)^(1-nu); %PMc
PMk=nu*(Ak)^(1-nu); %PMk
PK=Muk/Muc*Nk^(1-Muk)*PI/(Nc^(1-Muc)); %PK
PI=(Ak/Ac)^((1-nu)*gamma); %PI
