%% Steady-states
% Aqui são estados estacionários para diferentes pontos de operação
% uss = [50Hz, 50%];
% xss = [8311024.82175957;2990109.06207437;0.0099504224];
%----------------------------------------------------------
% uss = [35Hz, 100%]
% yss = [7466064.06187811 Pa, 289.27057711 m];
% xss = [9572888.40904742, 2123302.77208337, 0.00702289];
%----------------------------------------------------------
% uss = [35Hz, 10%]
% yss = [9258024.607735397 Pa, 326.379353557 m];
% xss = [11186036.43743129, 4690248.601938976, 0.00328039];

%% Simulation
close all
clear
clc
bcs_settings_hlim

% Importando dados
Entradas = readmatrix('Entradas_BCS_CSV.csv');

% Condição inicial 
x0 = [7536344.5283922757953405380249023, 3561883.2242104276083409786224365, 0.011747680694129930545344109305006]

% uk_1 = [0.001 0.001]';
uk_1 = [58 47]';    

Ts = 1; % período de amostragem
tsim = 1000; % tempo de simulação
nsim = tsim/Ts; % numero de amostras da simulação

for j=1:nsim
    
    tsim = j*Ts;
    uk_1 = Entradas(j,:);   
    [t,x_t]=ode45(@(t,x)bcs_model(t,x,uk_1) ,[0 Ts] ,x0); 
        % @(t,x) bcs_model(t,x,uk_1) é a função que descreve um sistema de EDO.
        % [0 Ts] intervalo de tempo
        % x0 valor inicial
        % uk_1[1] = fq; uk_1[2] = zc
        % x[1] = pbh ; x[2] = pwh; x[3] = qp
    xpk(j,:) = x_t(end,:);
    x0 = xpk(j,:);
    y_sea = eq_medicao(x0,uk_1);
    ypk = y_sea(1:2);
    pin(j,:) = ypk(1);
    H(j,:) = ypk(2);
    uk(j,:) = uk_1;
    Xk(:,j) = x0;
    % Potência
    Pk(j,:) = y_sea(end);
    
    
end
Saidas = [abs(pin), abs(H)];
%writematrix(Saidas, 'Saidas_xlsx.xlsx');
% Titulo
Colunas = {'P_in', 'H'};
Saidas = [Colunas; num2cell(Saidas)];
% Salva no Excel
writecell(Saidas, 'Saidas_xlsx.xlsx');

%% Simulation
close all
clear
clc
bcs_settings_hlim

% Importando dados
Entradas = readmatrix('Entradas_BCS_CSV.csv');

% Condição inicial 
x0 = [7536344.5283922757953405380249023, 3561883.2242104276083409786224365, 0.011747680694129930545344109305006]

% uk_1 = [0.001 0.001]';
uk_1 = [58 47]';    

Ts = 1; % período de amostragem
tsim = 1000; % tempo de simulação
nsim = tsim/Ts; % numero de amostras da simulação

for j=1:nsim
    
    tsim = j*Ts;
    uk_1 = Entradas(j,:);   
    [t,x_t]=ode45(@(t,x)bcs_model(t,x,uk_1) ,[0 Ts] ,x0); 
        % @(t,x) bcs_model(t,x,uk_1) é a função que descreve um sistema de EDO.
        % [0 Ts] intervalo de tempo
        % x0 valor inicial
        % uk_1[1] = fq; uk_1[2] = zc
        % x[1] = pbh ; x[2] = pwh; x[3] = qp
    xpk(j,:) = x_t(end,:);
    x0 = xpk(j,:);
    y_sea = eq_medicao(x0,uk_1);
    ypk = y_sea(1:2);
    pin(j,:) = ypk(1);
    H(j,:) = ypk(2);
    uk(j,:) = uk_1;
    Xk(:,j) = x0;
    % Potência
    Pk(j,:) = y_sea(end);
    
    
end
Saidas = [abs(pin), abs(H)];
%writematrix(Saidas, 'Saidas_xlsx.xlsx');
% Titulo
Colunas = {'P_in', 'H'};
Saidas = [Colunas; num2cell(Saidas)];
% Salva no Excel
writecell(Saidas, 'Saidas_xlsx.xlsx');
%% Grafico
t = Ts:Ts:nsim*Ts;
figure(1)
label = {'p_{bh} (Pa)','p_{wh} (Pa)','q_{p} (m^3/s)'};
for iy = 1:3
    subplot(3,1,iy)
    hold on
    plot(t,xpk(:,iy),'LineWidth',1.5) 
    ylabel(label(iy))
end
xlabel('Time (s)')

figure(2)
subplot(2,1,1)
hold on
grid on
plot(t,pin(:,1)/1e5,'LineWidth',1.5)
ylabel('p_{in} (bar)');
subplot(2,1,2)
hold on
grid on
plot(t,H(:,1),'LineWidth',1.5)
ylabel('H (m)');
xlabel('Time (s)')

label = {'f (Hz)','z_c (%)'};
figure(3)
for iu = 1:2
    subplot(2,1,iu)
    hold on
    grid on
    plot(t,uk(:,iu),'LineWidth',1.5)
    ylabel(label(iu))
end
xlabel('Time (s)')

% >>>>>>>>>>>>>>>>> Envelope colorbar <<<<<<<<<<<<<<<<<<<<
figure()
grid on
hold on
xx = Xk(3,:)*3600;
yy = H';
zz = zeros(size(Xk(3,:)));
time = Ts:Ts:nsim*Ts;
surface([xx;xx],[yy;yy],[zz;zz],[time;time],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',3);
    colormap jet
a = colorbar;
a.Label.String = 'Time (s)';
plot(Xk(3,1)*3600,H(1),'o','MarkerFaceColor',[0,1,0],'MarkerEdgeColor',[0,0,0])
plot(Xk(3,end)*3600,H(end),'o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0])
text(Xk(3,1)*3600,H(1),'t = 0','HorizontalAlignment','left')
text(Xk(3,end)*3600,H(end),sprintf('t = %d', tsim),'HorizontalAlignment','left')
BCS.Envelope.fig();
xlabel('q_p (m^3/h)')
ylabel('H (m)')

% Constantes
f0 = 60;
mu  = 0.025;  % Viscosity [Pa*s]
q0_dt = 25/3600; % Downtrhust flow at f0
q0_ut = 50/3600; % Uptrhust flow at f0
CH = -0.03*mu + 1;
Cq = 2.7944*mu^4 - 6.8104*mu^3 + 6.0032*mu^2 - 2.6266*mu + 1;
Cp = -4.4376*mu^4 + 11.091*mu^3 -9.9306*mu^2 + 3.9042*mu + 1;

%% Região de operação Downtrhust e Upthrust
H0_dt = -1.2454e6*q0_dt.^2 + 7.4959e3*q0_dt + 9.5970e2;
H0_dt = CH*H0_dt*(f0/f0).^2;
H0_ut = -1.2454e6*q0_ut.^2 + 7.4959e3*q0_ut + 9.5970e2;
H0_ut = CH*H0_ut*(f0/f0).^2;
% Variacao frequencia
f = linspace(30,70,1000); % Hz
H_ut = H0_ut*(f./f0).^2;
H_dt = H0_dt*(f./f0).^2;
% corrige lei da afinidade
Qdt = q0_dt.*f/f0;
Qut = q0_ut.*f/f0;
% Variacao vazao
flim = 35:5:65;
qop = linspace(0,q0_ut*flim(end)/f0,1000); % m3/s
Hop = zeros(length(flim),length(qop));
for i = 1:length(flim)
    q0 = qop./Cq*(f0/flim(i));
    H0 = -1.2454e6*q0.^2 + 7.4959e3*q0 + 9.5970e2;
    Hop(i,:) = CH*H0*(flim(i)/f0).^2;
end
% Calculo dos pontos de interseção para delimitação da região
[ip(1,1),ip(1,2)] = polyxpoly(qop*3600,Hop(1,:),Qdt*3600,H_dt);
[ip(2,1),ip(2,2)] = polyxpoly(Qdt*3600,H_dt,qop*3600,Hop(end,:));
[ip(3,1),ip(3,2)] = polyxpoly(qop*3600,Hop(end,:),Qut*3600,H_ut);
[ip(4,1),ip(4,2)] = polyxpoly(Qut*3600,H_ut,qop*3600,Hop(1,:));

% Ajuste do polinomio de frequencia maxima 65 Hz
p_35hz = polyfit(qop*3600,Hop(1,:),3);
H_35hz = @(qk) p_35hz*[cumprod(repmat(qk,length(p_35hz)-1,1),1,'reverse');ones(1,length(qk))];
q_35hz = linspace(ip(1,1),ip(4,1),100);
% Ajuste do polinomio de frequencia minima 35 Hz
p_65hz = polyfit(qop*3600,Hop(end,:),3);
H_65hz = @(qk) p_65hz*[cumprod(repmat(qk,length(p_65hz)-1,1),1,'reverse');ones(1,length(qk))];
q_65hz = linspace(ip(2,1),ip(3,1),100);
% Ajuste do polinomio de Downtrhust
p_dt = polyfit(Qdt*3600,H_dt,2);
H_dt = @(qk) p_dt*[cumprod(repmat(qk,length(p_dt)-1,1),1,'reverse');ones(1,length(qk))];
q_dt = linspace(ip(1,1),ip(2,1),100);
% Ajuste do polinomio de Uptrhust
p_ut = polyfit(Qut*3600,H_ut,2);
H_ut = @(qk) p_ut*[cumprod(repmat(qk,length(p_ut)-1,1),1,'reverse');ones(1,length(qk))];
q_ut = linspace(ip(4,1),ip(3,1),100);
% Constução da figura
BCS.Envelope.fig = @(aux) plot(q_35hz,H_35hz(q_35hz),':r',q_65hz,H_65hz(q_65hz),':r',q_ut,H_ut(q_ut),':r',q_dt,H_dt(q_dt),':r','LineWidth',2);
BCS.Envelope.ip = ip;
BCS.Envelope.fBounds = struct('H_35hz',H_35hz,'H_65hz',H_65hz,'H_dt',H_dt,'H_ut',H_ut);
% Função para a avaliação dos limites dada uma vazão.
BCS.Envelope.Hlim = @(qk) BoundHead(qk*3600,ip,BCS.Envelope.fBounds);

%% Subrotina
function Hlim = BoundHead(qk,ip,bounds)
    if qk < ip(1,1)
        Hlim = [ip(1,2),ip(1,2)];
    elseif qk < ip(2,1)
        Hlim = [bounds.H_35hz(qk);bounds.H_dt(qk)];
    elseif qk < ip(4,1)
        Hlim = [bounds.H_35hz(qk);bounds.H_65hz(qk)];
    elseif qk < ip(3,1)
        Hlim = [bounds.H_ut(qk);bounds.H_65hz(qk)];
    else
        Hlim = [ip(3,2),ip(3,2)];
    end
end
