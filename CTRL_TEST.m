% =================================================
% CONTROLADORES PARA BITBOT
% José Eduardo Molina Calzia 
% Carné 16063
% Agosto 2020
% =================================================


% Borrar estado anterior 
close all; 
clear all; 
% ================================================

% Selector de trayectoria o meta fija --> fg, traj
sel = 'fg';
%Opciones: pid, pidexp, nonlin, lqr, lqi
ctrl = 'nonlin'; 

graf = 1; 


%% Definición del sistema dinámico 
% xi = [x; y; theta]       u = [v;omega]


% Parámetros 
ell =0.01; %m

% Campo vectorial del sisrema dinámico 

f = @(xi,u)[u(1)*cos(xi(3)); ...
            u(1)*sin(xi(3)); ...
            u(2)]; 
        
        
% Difeomorfismo
finv = @(xi,mu) [1,0; 0,1/ell] * [cos(xi(3)), -sin(xi(3)); sin(xi(3)), cos(xi(3))]' * [mu(1); mu(2)];
% Matrices del sistema linealizado 

A = zeros(2); 
B = eye(2); 


%% Parámetros de simulación 
dt = 0.001; 
t0 = 0; 
tf = 15; 
N = (tf-t0)/dt; 

 

%% Condiciones iniciales

%xi0 = zeros(3,1); 
xi0 = [-4;-5;0]; 
u0 = zeros(2,1); 

% Vector de estado
xi = xi0; 
% Vector de entradas
u = u0; 

% Array de históticos
XI = zeros(numel(xi),N+1);
U = zeros(numel(u),N+1);
XI(:,1) = xi; 
U(:,1) = u; 

%% Meta fija


xg = 5; 
yg = 5; 
thetag = pi; 


%% Trayectoria
 
t = 0:dt:tf;
traj = 2*[cos(0.1*2*pi*t); sin(0.1*2*pi*t)];


%% Controladores 




% PID 
    % Posición
kpP = 1; 
kiP = 0.001; 
kdP = 0.1; 
EP = 0; 
eP_1 = 0; 

    % Orientación 
kpO = 1;
kiO = 0.001; 
kdO = 0.1;
EO = 0;
eO_1 = 0;

    % Acercamiento exp 
v0 = 2; 
alpha = 0.5; 

% Controlador no lineal 
k_rho = 10; 
k_alpha = 20; 
k_beta = -10; 

% LQR
Q = eye(2); 
Q(1,1) = 0.0001; 
Q(2,2) = 0.0001; 

R = eye(2); 
R(1,1) = 1; 
R(2,2) = 1; 

    % Calculo de Klqr
Klqr = lqr(A,B,Q,R); 


% LQI
Cr = eye(2); 
Dr = zeros(2); 

AA = [A, zeros(size(Cr')); Cr, zeros(size(Cr,1))];
BB = [B; Dr];

QQ = eye(size(A,1) + size(Cr,1));
QQ(3,3) = 1; 
QQ(4,4) = 1; 

RR = eye(2); 
RR(1,1) = 1; 
RR(2,2) = 1; 

    %Calculo de Klqr
Klqi = lqr(AA,BB,QQ,RR)
ref = [xg;yg]; 
sigma = 0; 



%% Solución del sistema dinámico

for n = 0:N
    if (strcmp(sel,'traj'))
        xg = traj(1,n+1);
        yg = traj(2,n+1);
    end
    
    
    if(strcmp(ctrl, 'pid'))
        fprintf("PID sin exp\n\n"); 
        x = xi(1); 
        y = xi(2); 
        theta = xi(3);
        
        % Error
        e = [xg-x; yg-y];
        thetag = atan2(e(2),e(1)); 
        eP = norm(e); 
        eO = thetag-theta; 
        eO = atan2(sin(eO), cos(eO));
        
        % Pos --> Velocidad lineal
        eP_D = eP-eP_1;
        EP = EP + eP;
        v = kpP*eP + kiP*EP + kdP*eP_D;
        eP_1 = eP;
        % Ori --> Velocidad angular 
        eO_D = eO - eO_1;
        EO = EO + eO;
        w = kpO*eO + kiO*EO + kdO*eO_D;
        eO_1 = eO;      
        
        %Vector de entradas
        u = [v; w];

    elseif(strcmp(ctrl, 'pidexp'))
        fprintf("PID con exp\n\n"); 
        x = xi(1); 
        y = xi(2); 
        theta = xi(3);
        
       % Error 
       e = [xg-x; yg-y];
       thetag = atan2(e(2),e(1)); 
       eP = norm(e);
       eO = thetag - theta;
       eO = atan2(sin(eO), cos(eO));
      
       % Pos --> Velocidad lineal
       constP = v0 * (1-exp(-alpha*eP^2)) / eP;
       v = constP*eP;
       
       % Ori --> Velocidad angular 
       eO_D = eO - eO_1;
       EO = EO + eO;
       w = kpO*eO + kiO*EO + kdO*eO_D;
       eO_1 = eO;       
       
       % Vector de entradas
       u = [v; w];

  
    elseif(strcmp(ctrl, 'nonlin'))
       fprintf("No linear\n\n"); 
       x = xi(1); 
       y = xi(2); 
       theta = xi(3);
       e = [xg - x; yg - y];
            
       rho = norm(e);
       alpha = -theta + atan2(e(2), e(1));
       beta = -theta - alpha;
            
       v = k_rho*rho;
       w = k_alpha*alpha + k_beta*beta;
            
       % Vector de entradas
       u = [v;w]; 
        

    elseif(strcmp(ctrl, 'lqr'))
       fprintf("LQR\n\n"); 
       x = xi(1); 
       y = xi(2);
       e = [x - xg; y - yg];
       mu = -Klqr*e;

       % Vector de entradas 
       u = finv(xi, mu);
    elseif(strcmp(ctrl, 'lqi'))
        fprintf("LQI\n\n"); 
        x = xi(1:2); 
        ref = [xg;yg]; 
        sigma = sigma + (Cr*x - ref)*dt;
        mu = -Klqi*[x; sigma];
        
        % Vector de entradas
        u = finv(xi, mu); 
        
    else %default pid 
        fprintf("DEF\n\n"); 
        x = xi(1); 
        y = xi(2); 
        theta = xi(3);
        
        % Error
        e = [xg-x; yg-y];
        thetag = atan2(e(2),e(1)); 
        eP = norm(e); 
        eO = thetag-theta; 
        eO = atan2(sin(eO), cos(eO));
        
        % Pos --> Velocidad lineal
        eP_D = eP-eP_1;
        EP = EP + eP;
        v = kpP*eP + kiP*EP + kdP*eP_D;
        eP_1 = eP;
        % Ori --> Velocidad angular 
        eO_D = eO - eO_1;
        EO = EO + eO;
        w = kpO*eO + kiO*EO + kdO*eO_D;
        eO_1 = eO;      
        
        %Vector de entradas
        u = [v; w];
    end
    
    % RK4 --> actualización del sistema dinámico
    k1 = f(xi, u);
    k2 = f(xi+(dt/2)*k1, u);
    k3 = f(xi+(dt/2)*k2, u);
    k4 = f(xi+dt*k3, u);
    xi = xi + (dt/6)*(k1+2*k2+2*k3+k4);
    
    
    % Se guardan las trayectorias del estado y las entradas
    XI(:,n+1) = xi;
    U(:,n+1) = u;
    
    
end

if(graf == 1)
%% Gráficas
figure('units','normalized','outerposition',[0,0,1,1]);
t = t0:dt:tf; 
plot(t,XI'); 
xlabel('$t[s]$', 'Interpreter', 'latex', 'Fontsize', 15); 
ylabel('$x(t)$', 'Interpreter', 'latex', 'Fontsize', 15); 
l = legend('$x(t)$', '$y(t)$', '$\theta(t)$');
set(l, 'Interpreter', 'latex', 'FontSize', 12);
grid minor
if(strcmp(ctrl, 'pid'))
    title('PID','Interpreter', 'latex', 'Fontsize', 25 )
elseif(strcmp(ctrl, 'pidexp'))
    title('PID acercamiento exponencial','Interpreter', 'latex', 'Fontsize', 25 )
elseif(strcmp(ctrl, 'lqi'))
    title('LQI','Interpreter', 'latex', 'Fontsize', 25 )
elseif(strcmp(ctrl, 'lqr'))
    title('LQR','Interpreter', 'latex', 'Fontsize', 25 )
elseif(strcmp(ctrl, 'nonlin'))
    title('No linear','Interpreter', 'latex', 'Fontsize', 25 )
end


figure('units','normalized','outerposition',[0,0,1,1]);
t = t0:dt:tf; 
plot(t,U'); 
xlabel('$t[s]$', 'Interpreter', 'latex', 'Fontsize', 15); 
ylabel('$u(t)$', 'Interpreter', 'latex', 'Fontsize', 15); 
l = legend('$v(t)$', '$w(t)$');
set(l, 'Interpreter', 'latex', 'FontSize', 12);
grid minor
if(strcmp(ctrl, 'pid'))
    title('PID','Interpreter', 'latex', 'Fontsize', 25 )
elseif(strcmp(ctrl, 'pidexp'))
    title('PID acercamiento exponencial','Interpreter', 'latex', 'Fontsize', 25 )
elseif(strcmp(ctrl, 'lqi'))
    title('LQI','Interpreter', 'latex', 'Fontsize', 25 )
elseif(strcmp(ctrl, 'lqr'))
    title('LQR','Interpreter', 'latex', 'Fontsize', 25 )
elseif(strcmp(ctrl, 'nonlin'))
    title('No linear','Interpreter', 'latex', 'Fontsize', 25 )
end





%% Simulación visual
figure('units','normalized','outerposition',[0,0,1,1]);
s = max(max(abs(XI(1:2,:))));
xlim(s*[-1, 1]+[-2, 2]);
ylim(s*[-1, 1]+[-2, 2]);
grid minor;
hold on;

q = XI(:,1);
x = q(1); 
y = q(2); 
theta = q(3);

if(strcmp(sel,'traj'))
    plot(traj(1,:), traj(2,:),'k'); 
elseif (strcmp(sel,'fg'))
    plot(xg,yg,'o'); 
end


trajplot = plot(x, y, '--k', 'LineWidth', 1, 'Color', [0,1,0]);
p1 = [0-3.5;4];p2 = [3-3.5;4];p3 = [3-3.5;2];p4 = [7-3.5;2];
p5 = [7-3.5;-2];p6 = [3-3.5;-2];p7 = [3-3.5;-4];p8 = [0-3.5;-4];
BBt =[p1,p2,p3,p4,p5,p6,p7,p8]/20;
bodyplot = fill(BBt(1,:)+x,BBt(2,:)+y,[1,0,0]); 

xlabel('$x$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$y$', 'Interpreter', 'latex', 'Fontsize', 16);

if(strcmp(ctrl, 'pid'))
    title('PID','Interpreter', 'latex', 'Fontsize', 25 )
elseif(strcmp(ctrl, 'pidexp'))
    title('PID acercamiento exponencial','Interpreter', 'latex', 'Fontsize', 25 )
elseif(strcmp(ctrl, 'lqi'))
    title('LQI','Interpreter', 'latex', 'Fontsize', 25 )
elseif(strcmp(ctrl, 'lqr'))
    title('LQR','Interpreter', 'latex', 'Fontsize', 25 )
elseif(strcmp(ctrl, 'nonlin'))
    title('No linear','Interpreter', 'latex', 'Fontsize', 25 )
end


hold off;



for n = 2:N+1
    q = XI(:,n);
    x = q(1); y = q(2); theta = q(3);
    
    trajplot.XData = [trajplot.XData, x];
    trajplot.YData = [trajplot.YData, y];
    
    rot = [cos(theta-0), -sin(theta-0); sin(theta-0), cos(theta-0)] * BBt;
    bodyplot.XData = rot(1,:) + x;
    bodyplot.YData = rot(2,:) + y;
  
    pause(dt/2);
end
end
