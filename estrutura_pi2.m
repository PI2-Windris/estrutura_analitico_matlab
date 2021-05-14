%% ALGORITMO PARA SOLUÇÃO ESTRUTURAL
%%% GRUPO 3 - ESTAÇÃO WINDRIS
%%% Universidade de Brasília
%%% Faculdade do Gama

%% Parte 1 - Definição do aerofólio e cálculo da eficiência
clc, clear all

V=5;    %velocidade do vento [m/s]
TSR=2.64;   %tip speed ratio
sigma=0.24; %solidez
n=3;        %número de pás
rho=1.23;   %densidade do ar [kg/m³]
AR=7.14;    %razão de aspecto

%Dados geométricos da turbina
R=0.6;  %raio do rotor [m]
d=2*R;  %diâmetro do rotor [m]
b=2.4*R;    %comprimento das pás [m]
c=b/AR;     %comprimento de corda do aerofólio [m]

%velocidade tangencial [m/s]
Vt=TSR*V;   %Vt=Rw --> w=Vt/R --- TSR*V=Rw

%Theta e alfa em graus
theta=[1:1:360];    %ângulo azimutal [deg]
alfa=atand(sind(theta)./(TSR+cosd(theta))); %ângulo de ataque [deg]

%Ângulos em radianos
alfar=alfa*pi/180;  %rad
beta=asin(Vt/V.*sin(alfar));    %rad 
gama=(180-alfar-beta);          %rad

%Cálculo da velocidade relativa
W=sqrt(V^2+Vt^2-2*V*Vt*cos(gama));  %[m/s]

%Coeficientes aerodinamicos da simulação numérica (caso de TSR=2.64)
naca0015=xlsread('naca0015_TSR2');      selig1210=xlsread('selig1210_TSR2');    %leitura das tabelas com os resultados

cl15=naca0015(:,2);     cd15=naca0015(:,3);     cm15=naca0015(:,4); %coeficiente de sustentação, arrasto e momento para o NACA 0015
cl12=selig1210(2:361,2);     cd12=selig1210(2:361,3);     cm12=selig1210(2:361,4); %coeficiente de sustentação, arrasto e momento para o Selig 1210

%pressão dinâmica [kg/m²s]
q=0.5*rho*W.^2;  

%área do aerofólio [m²]
S=b*c;

for i=1:length(theta)
    %forças aerodinamicas para o NACA 0015
    L15(i)=q(i)*S*cl15(i);
    D15(i)=q(i)*S*cd15(i);
    N15(i)=L15(i)*cosd(alfa(i))+D15(i)*sind(alfa(i));
    A15(i)=-L15(i)*sind(alfa(i))+D15(i)*cosd(alfa(i));

    %Calculo do torque para o NACA 0015
    T15(i)=R*[N15(i)*(sind(theta(i)))+A15(i)*(cosd(theta(i)))]; %N*m
    
    %forças aerodinamicas para o Selig 1210
    L12(i)=q(i)*S*cl12(i);
    D12(i)=q(i)*S*cd12(i);
    N12(i)=L12(i)*cosd(alfa(i))+D12(i)*sind(alfa(i));
    A12(i)=-L12(i)*sind(alfa(i))+D12(i)*cosd(alfa(i));

    %Calculo do torque para o Selig 1210
    T12(i)=R*[N12(i)*(sind(theta(i)))+A12(i)*(cosd(theta(i)))]; %N*m
end

ct15=(-1)*T15./(q*S*R); %coeficiente de torque médio para o NACA 0015
Tmedio15=sum(ct15*R*0.5*rho*S*V^2)/(2*pi);  %torque médio para o NACA 0015
Pmedio15=Tmedio15*Vt/R; %potência média para o NACA 0015

ct12=T12./(q*S*R); %coeficiente de torque médio para o Selig 1210
Tmedio12=sum(ct12*R*0.5*rho*S*V^2)/(2*pi);  %torque médio para o Selig 1210
Pmedio12=Tmedio12*Vt/R; %potência média para o Selig 1210

%===================== Resultados ====================%
tabela1=table([theta'],[alfa'],[W'],'VariableNames',{'theta' 'alfa' 'Velocidade_Relativa'});
tabela2=table([Tmedio15],[Tmedio12],[Pmedio15],[Pmedio12],'VariableNames',{'TmedioNACA' 'TmedioSelig' 'PmediaNACA' 'PmediaSelig'})

%===================== Gráficos ====================%
figure(1)
plot(theta,T15,'r',theta,T12,'b')
grid on
legend('NACA 0015','Seling S1210')
xlabel('Ângulo azimutal [graus]')
ylabel('Torque instantâneo [Nm]')

figure(2)
plot(theta,alfa,'b')
grid on
xlabel('Ângulo azimutal [graus]')
ylabel('Ângulo de ataque [graus]')

figure(3)
plot(theta,cl15,'r',theta,cl12,'b')
grid on
legend('NACA 0015','Seling S1210')
xlabel('Ângulo azimutal [graus]')
ylabel('c_l')

figure(4)
polarplot(cd15,'LineWidth',2)
hold on
polarplot(cd12,'LineWidth',2)
legend('NACA 0015','Seling S1210');
ax=gca;
ax.ThetaZeroLocation = 'top';
ax.RTick = [0:0.5:2];
ax.ThetaMinorGrid = 'on';
ax.RMinorTick = 'on';
hold off

figure(5)
plot(theta,cd15,theta,cd12)
legend('NACA0015','Selig1210')
xlabel('Ângulo azimutal')
ylabel('c_d')

figure(6)
plot(theta,cm15,theta,cm12)
legend('NACA0015','Selig1210')
xlabel('Ângulo azimutal')
ylabel('c_m')

%% Parte 2 - Análise da variação da potência/torque em função da velocidade do vento

TSR=2.64;   %tip speed ratio
sigma=0.24; %solidez
n=3;        %número de pás
rho=1.23;   %densidade do ar [kg/m³]
AR=7.14;    %razão de aspecto

%Dados geométricos da turbina
R=0.6;  %raio do rotor [m]
d=2*R;  %diâmetro do rotor [m]
b=2.4*R;    %comprimento das pás [m]
c=b/AR;     %comprimento de corda do aerofólio [m]

for V=1:6
    %velocidade tangencial [m/s]
    Vt=TSR*V;
    
    %velocidade angular [rad/s]
    w=Vt/R;

    %Theta e alfa em graus
    theta=[1:1:360];    %ângulo azimutal [deg]
    alfa=atand(sind(theta)./(TSR+cosd(theta))); %ângulo de ataque [deg]

    %Ângulos em radianos
    alfar=alfa*pi/180;  %rad
    beta=asin(Vt/V.*sin(alfar));    %rad 
    gama=(180-alfar-beta);          %rad

    %Cálculo da velocidade relativa
    W=sqrt(V^2+Vt^2-2*V*Vt*cos(gama));  %[m/s]

    %Coeficientes aerodinamicos da simulação numérica (caso de TSR=2.64)
    naca0015=xlsread('naca0015_TSR2');      selig1210=xlsread('selig1210_TSR2');    %leitura das tabelas com os resultados

    cl15=naca0015(:,2);     cd15=naca0015(:,3);     cm15=naca0015(:,4); %coeficiente de sustentação, arrasto e momento para o NACA 0015
    cl12=selig1210(2:361,2);     cd12=selig1210(2:361,3);     cm12=selig1210(2:361,4); %coeficiente de sustentação, arrasto e momento para o Selig 1210

    %pressão dinâmica [kg/m²s]
    q=0.5*rho*W.^2;  

    %área do aerofólio [m²]
    S=b*c;

    for i=1:length(theta)
        %forças aerodinamicas para o NACA 0015
        L15(i)=q(i)*S*cl15(i);
        D15(i)=q(i)*S*cd15(i);
        N15(i)=L15(i)*cosd(alfa(i))+D15(i)*sind(alfa(i));
        A15(i)=-L15(i)*sind(alfa(i))+D15(i)*cosd(alfa(i));

        %Calculo do torque para o NACA 0015
        T15(i)=R*[N15(i)*(sind(theta(i)))+A15(i)*(cosd(theta(i)))]; %N*m

        %forças aerodinamicas para o Selig 1210
        L12(i)=q(i)*S*cl12(i);
        D12(i)=q(i)*S*cd12(i);
        N12(i)=L12(i)*cosd(alfa(i))+D12(i)*sind(alfa(i));
        A12(i)=-L12(i)*sind(alfa(i))+D12(i)*cosd(alfa(i));

        %Calculo do torque para o Selig 1210
        T12(i)=R*[N12(i)*(sind(theta(i)))+A12(i)*(cosd(theta(i)))]; %N*m
    end
    ct15=(-1)*T15./(q*S*R); %coeficiente de torque médio para o NACA 0015
    Tmedio15(V)=sum(ct15*R*0.5*rho*S*V^2)/(2*pi); %torque médio para o NACA 0015
    Pmedio15(V)=Tmedio15(V)*w; %potência média para o NACA 0015

    ct12=T12./(q*S*R); %coeficiente de torque médio para o Selig 1210
    Tmedio12(V)=sum(ct12*R*0.5*rho*S*V^2)/(2*pi); %torque médio para o Selig 1210
    Pmedio12(V)=Tmedio12(V)*w; %potência média para o Selig 1210
    
    vel(V)=V;   %vetor para construção dos gráficos
end

figure(7)
plot(vel,Tmedio12,'r')
grid on
xlim([1 6])
xlabel('Velocidade do vento [m/s]')
ylabel('Torque médio [Nm]')

figure(8)
plot(vel,Pmedio12,'r')
grid on
xlim([1 6])
xlabel('Velocidade do vento [m/s]')
ylabel('Potência média [W]')