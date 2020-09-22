% MATLAB controller for Webots
% File:          Kalman.m
% Date:
% Description:
% Author:
% Modifications:

% uncomment the next two lines if you want to use
% MATLAB's desktop to interact with the controller:
desktop;
%keyboard;

clear all; 
close all; 



TIME_STEP = 50;
WHEEL_RADIUS = 195/2000; 
DISTANCE_CENTER = 381/2000;
MAX_SPEED = 12.3; 
left_wheel = wb_robot_get_device('left wheel');
right_wheel = wb_robot_get_device('right wheel');
position_sensor = wb_robot_get_device('gps');
orientation_sensor = wb_robot_get_device('compass');

wb_gps_enable(position_sensor, TIME_STEP);
wb_compass_enable(orientation_sensor, TIME_STEP);
encoder_R = wb_robot_get_device('right wheel sensor'); 
encoder_L = wb_robot_get_device('left wheel sensor'); 

wb_motor_set_position(left_wheel, inf);
wb_motor_set_position(right_wheel, inf);
wb_motor_set_velocity(left_wheel, 0.0);
wb_motor_set_velocity(right_wheel, 0.0);
wb_gps_enable(position_sensor, TIME_STEP);
wb_compass_enable(orientation_sensor, TIME_STEP);

wb_position_sensor_enable(encoder_R, TIME_STEP)
wb_position_sensor_enable(encoder_L, TIME_STEP)





robot_r =  (195/1000)/2; % radio de las ruedas
robot_ell = 381/1000; % distancia entre ruedas
robot_N = 256; % pulses por revolución de los encoders




desvx = 0.07/100; 
desvy = 0.06/100; 
vartheta = 1.03/10; %grados  
desvtheta = vartheta*pi/180; % radianes






%==================================================
%&%%%%% INICIALIZACIÓN 

% Condiciones iniciales 
delta0 = 0; 

x0 = 2 + delta0; 
y0 = 0 + delta0; 
theta0 = 0; 


xi0 = [x0; y0; theta0];
mu0 = [0; 0];
xi = xi0; 
mu = mu0; 

% Arrays para almacenar trayectorias 

XI = xi0; 
MU = mu0; 


% Arrays auxiliares 
CAM = [0;0;0]; 

ENCR = [0]; 
ENCL = [0]; 


%============================================================================
 
%%%%% DEFINICIÓN DE SISTEMAS DINÁMICOS  

% Modelo de odometría 
odo = @(x,drho,dtheta) [x(1)+drho*cos(x(3)); x(2)+drho*sin(x(3)); x(3)+dtheta]; 


% Jacobianos del modelo de odometría
dodo_dx = @(x,drho,dtheta) [1,0, -drho*sin(x(3)); 0,1, drho*cos(x(3)); 0,0,1];
dodo_dw = @(x,drho,dtheta) [cos(x(3)),0; sin(x(3)),0; 0,1];


Cd = [1,0,0; 0,1,0; 0,0,1]; 



% Varianzas del encoder (proceso)
max_speed = 1.2; %m/s
delta_vmax = max_speed*TIME_STEP/1000; 
sw1 = delta_vmax; 
sw3 = 2*pi/180; % se sobre-estima el error angular
Qw = diag([sw1^2, sw3^2]); 


% Varianzas del GPS (observación)
sv1 = desvx; 
sv2 = desvy; 
sv3 = desvtheta; 
Qv = diag([sv1^2, sv2^2, sv3^2]); 


% Mediciones previas 
pulsesR = 0; 
pulsesL = 0; 

DR = 0; 
DL = 0; 

drho = 0; %pos lineal

% Estimados a-prior y a-post
xhat0 = [x0; y0; theta0];
sys_order = numel(xhat0);
xhat_prior = zeros(3,1); 
xhat_post = xhat0;


% Matriz de coviarianza
sigma_e = 0.001; 
P_prior = zeros(3, 3); 
P_post =eye(3) * sigma_e^2; 


XHAT = xhat0; 




cam_count = 0; 
y_1 = 0; 



%=============================================================
% Controlador 

xg = 7; 
yg = 5;
thetag = pi; 

ctrl = 'lqr'; 
sel = 1; 


% Parámetros de simulación 


if(strcmp(ctrl, 'pid'))
tf = 25; 
elseif(strcmp(ctrl, 'pidexp'))
tf = 25; 
elseif(strcmp(ctrl, 'lqi'))
tf = 200;  
elseif(strcmp(ctrl, 'lqr'))
tf = 100; 
elseif(strcmp(ctrl, 'nonlin'))
tf = 12.5; 
end

if(sel)
  tf = 1000; 
end
t0 = 0; 
%tf = 12.5; 
count = 1; 
tmax = (tf-t0)/(TIME_STEP/1000); 




%% Controlador no lineal 

k_rho = 5; 
k_alpha = 8; 
k_beta = -1.5; 

u0 = [0;0]; 
U = u0; 

% LQI y LQR

u1 = [0;0]; 
Ekx = 0; ekx = 0; 
Eky = 0; eky = 0; 
v = 0; 
w = 0; 
dt = TIME_STEP/1000;
ell = 0.1; 
difeo = @(theta,mu) [1,0;0,1/ell]*[cos(theta), -sin(theta); sin(theta), cos(theta)]'*[mu(1); mu(2)]; 




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

%% Acercamiento exp 
alpha = 0.5; 


%%%%%% TRAJ 
t = 0:dt:tf;
splinex = 2:0.25:7.5;  
spliney = [0, 0.5, 1.5, 2.3, 3, 3.5, 3.8, 4.1, 4.3, 4.6, 4.7, 4.75, 4.65, 4.3, 4.25, 4.15, 3.8, 3.3, 2.75, 2, 1.25, 0.5, 0]; 


traj = [splinex;spliney]; 
er = [0;0]; 
if(sel)
  xg = traj(1,1); 
  yg = traj(2,1); 
end
ct = 1; 


% Loop de simulación 
while wb_robot_step(TIME_STEP) ~= -1

  
  %============================================ 
  %%%%%% ESTIMACIÓN DE POSICIÓN Y ÁNGULO  
  
  % Mediciones y Varianzas de camara 
  ran = randn(3,1); 
  ranx = ran(1)*desvx; 
  rany = ran(2)*desvy; 
  rantheta = ran(3)*desvtheta; 
  
  cam_pos = wb_gps_get_values(position_sensor)'; 
  cam_x = cam_pos(1)+ranx; 
  cam_y = cam_pos(3)+rany; 
  cam_position= [cam_x; cam_y]; 
    
  norte = wb_compass_get_values(orientation_sensor);
  cam_theta_prim  = atan2(norte(1), norte(3));
  cam_theta = cam_theta_prim + rantheta; 
  
  
  % Mediciones y Varianzas de encoders 
  encr = wb_position_sensor_get_value(encoder_R); 
  encl = wb_position_sensor_get_value(encoder_L); 
  desv_pulses = randn(1,1)*(max_speed*robot_r*TIME_STEP/1000); 
  
  encr_pulses = encr*robot_N/(2*pi)+desv_pulses; 
  encl_pulses = encl*robot_N/(2*pi)+desv_pulses; 
  
  
  
  xi = [cam_pos(1); cam_pos(3); cam_theta_prim]; 
  
  %%%% EKF 
  % Incrementos de posición lineal y angular --> modelo de odometría 
  DL = 2*pi*robot_r*((encl_pulses - pulsesL)/robot_N);
  DR = 2*pi*robot_r*((encr_pulses - pulsesR)/robot_N);
  
  pulsesL = encl_pulses; 
  pulsesR = encl_pulses; 
  
  drho = (DL+DR)/2; 
  dtheta = (DR - DL)/robot_ell;
  
  
  % Jacobianos 
  Ad = dodo_dx(xhat_post, drho, dtheta);
  Fd = dodo_dw(xhat_post, drho, dtheta); 
  
  % Predicción
  xhat_prior = odo(xhat_post, drho, dtheta);
  P_prior = Ad*P_post*Ad' + Fd*Qw*Fd';
  
  % Corrección 
  
  %{
  y = [cam_position; cam_theta]; 
  CAM_SENSOR = y; 
  L_k = P_prior*Cd'*(Qv + Cd*P_prior*Cd')^(-1);
  xhat_post = xhat_prior + L_k*(y - Cd*xhat_prior);
  P_post = P_prior - L_k*Cd*P_prior;
  y_1 = y; 
  %}
  
  
  

  if(mod(cam_count,2) == 0)
  %fprintf("Corrección\n")
  y = [cam_position; cam_theta]; 
  CAM_SENSOR = y; 
  L_k = P_prior*Cd'*(Qv + Cd*P_prior*Cd')^(-1);
  xhat_post = xhat_prior + L_k*(y - Cd*xhat_prior);
  P_post = P_prior - L_k*Cd*P_prior;
  y_1 = y; 
  else
  y = [cam_position; cam_theta]; 
  CAM_SENSOR = y;
  y = y_1; 
  %xhat_post = xhat_prior; 
  %P_post = P_prior; 
  L_k = P_prior*Cd'*(Qv + Cd*P_prior*Cd')^(-1);
  xhat_post = xhat_prior + L_k*(y - Cd*xhat_prior);
  P_post = P_prior - L_k*Cd*P_prior;
  end 
  cam_count = cam_count+1; 
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CONTROLADOR DE POSE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(sel)
    er = [xg - xhat_post(1); yg - xhat_post(2)];
    if (norm(er)<0.05)
      ct = ct+1; 
      if(ct>size(traj,2))
      break
      end
      xg = traj(1,ct);
      yg = traj(2,ct);
      fprintf("Cambio\n"); 
    end
  end
  
  
  
  if(strcmp(ctrl, 'nonlin'))
      theta = xhat_post(3);
      e = [xg - xhat_post(1); yg - xhat_post(2)];
            
      rho = norm(e);
      alpha = -theta + atan2(e(2), e(1));
      alpha = atan(sin(alpha)/cos(alpha)); 
      beta = -theta - alpha;
            
      v = k_rho*rho;
      w = k_alpha*alpha + k_beta*beta;
      u = [v;w];
  
      speedR = (v-w*DISTANCE_CENTER/WHEEL_RADIUS);
      speedL = (v+w*DISTANCE_CENTER/WHEEL_RADIUS);
  end 
  
  if(strcmp(ctrl,'lqi'))
      theta = xhat_post(3); 
 
      x = xhat_post(1:2); 
      ref = [xg;yg]; 
        
        
      ekx = x(1) - ref(1); 
      Ekx = Ekx+ekx; 
        
      eky = x(2) - ref(2); 
      Eky = Eky+eky; 
        
        
      u1(1) = -10*sqrt(3)*x(1) - 1*Ekx*dt; 
      u1(2) = -10*sqrt(3)*x(2) - 1*Eky*dt; 
      
      u = -difeo(theta,u1);
      v = u(1); w = u(2); 
  
      speedL = (v-w*DISTANCE_CENTER/WHEEL_RADIUS)/50;
      speedR = (v+w*DISTANCE_CENTER/WHEEL_RADIUS)/50;

  end
  
  
  
   if(strcmp(ctrl, 'pid'))
        x = xhat_post(1); 
        y = xhat_post(2); 
        theta = xhat_post(3);
        
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
        speedL = (v-w*DISTANCE_CENTER/WHEEL_RADIUS);
        speedR = (v+w*DISTANCE_CENTER/WHEEL_RADIUS);
        
        
  elseif(strcmp(ctrl, 'pidexp'))
    
        x = xhat_post(1); 
        y = xhat_post(2); 
        theta = xhat_post(3);
        
       % Error 
       e = [xg-x; yg-y];
       thetag = atan2(e(2),e(1)); 
       eP = norm(e);
       eO = thetag - theta;
       eO = atan2(sin(eO), cos(eO));
      
       % Pos --> Velocidad lineal
       constP = MAX_SPEED * (1-exp(-alpha*eP^2)) / eP;
       v = constP*eP;
       
       % Ori --> Velocidad angular 
       eO_D = eO - eO_1;
       EO = EO + eO;
       w = kpO*eO + kiO*EO + kdO*eO_D;
       eO_1 = eO;       
       
       % Vector de entradas
       u = [v; w];
       speedL = (v-w*DISTANCE_CENTER/WHEEL_RADIUS);
       speedR = (v+w*DISTANCE_CENTER/WHEEL_RADIUS) ;
   elseif(strcmp(ctrl, 'lqr'))
  
       x = xhat_post(1); 
       y = xhat_post(2);
       e = [x - xg; y - yg];
         
       u1(1) = -5*e(1); 
       u1(2) = -5*e(2);  

       % Vector de entradas 
       u = -difeo(xhat_post(3),u1); 
       v = u(1); w = u(2); 
       speedL = (v-w*DISTANCE_CENTER/WHEEL_RADIUS)/10;
       speedR = (v+w*DISTANCE_CENTER/WHEEL_RADIUS)/10;
  end 
       
  
  wb_motor_set_velocity(left_wheel, speedL); 
  wb_motor_set_velocity(right_wheel, speedR);

  
   

  % Guardar salida del EKF
  XHAT = [XHAT, xhat_post]; 
  
  % Trayectorias de estado reales 
  XI = [XI,xi]; 
  MU = [MU,mu]; 
  
  % Mediciones de sensores
  
  CAM = [CAM, CAM_SENSOR]; 
  ENCR = [ENCR, encr_pulses]; 
  ENCL = [ENCL, encl_pulses]; 
  
  
  % Trayectorias del controlador 
  
  U = [U,u]; 
    
  
  %wb_console_print(sprintf('%2.2f x %2.2f\n\n', size(P_prior(1)), size(P_prior(2))),WB_STDOUT);

  % if your code plots some graphics, it needs to flushed like this:
  drawnow;
  
 
  count = count + 1; 
  if (count > tmax)
  break; 
  end
  

end
t = linspace(t0,tf,count); 
t2 = t0:TIME_STEP/1000:tf;

disp(size(t2)); 
disp(size(XI));  


figure(1); 
plot(t, XI', 'LineWidth', 1); 
hold on; 
ax = gca; 
%ax.ColorOrderIndex = 1; 
plot(t, XHAT', '--','LineWidth', 1); 
hold off; 
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\mathbf{x}(t)$', 'Interpreter', 'latex', 'Fontsize', 16);
yline(yg, '--', '$y_g$','Interpreter', 'latex'); 
yline(xg, '-.', '$x_g$','Interpreter', 'latex'); 
l = legend('$x(t)$', '$y(t)$', '$\theta(t)$','$\hat{x}(t)$','$\hat{y}(t)$','$\hat{\theta}(t)$',...
 'Location', 'best', ...
    'Orientation', 'vertical');
set(l, 'Interpreter', 'latex', 'FontSize', 12);
grid minor;
title('Respuesta del controlador', 'Interpreter', 'latex', 'Fontsize', 20)
ylim([-pi*1.2, 8.5])

figure(2);
plot(t, CAM', 'LineWidth', 1);
hold on;
ax = gca; 
%ax.ColorOrderIndex = 1; 
plot(t, XHAT','--' ,'LineWidth', 1);
hold off;
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\mathrm{pos}(t)$', 'Interpreter', 'latex', 'Fontsize', 16);
title('Comparacion salida EKF vs CAM');
l = legend('$x_{cam}(t)$', '$y_{cam}(t)$', '$\theta_{cam}(t)$','$x_{EKF}(t)$', '$y_{EKF}(t)$', '$\theta_{EKF}(t)$',...
 'Location', 'best','Orientation', 'vertical');
set(l, 'Interpreter', 'latex', 'FontSize', 12);
grid minor;


if(sel)
figure(3)
plot(traj(1,:),traj(2,:), 's', 'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6], 'MarkerSize', 10);
hold on 
plot(x0,y0,'o', 'MarkerEdgeColor','blue','MarkerFaceColor',[.6 .6 1], 'MarkerSize', 10);
hold on
plot(XI(1,:)',XI(2,:)','--k', 'LineWidth',1)
title('Trayectoria del robot', 'Interpreter', 'latex', 'Fontsize', 20); 
xlabel('$x$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$y$', 'Interpreter', 'latex', 'Fontsize', 16);
l = legend('$\zeta_g$','$\zeta_0$','Trayectoria','Location', 'best','Orientation', 'vertical');
set(l, 'Interpreter', 'latex', 'FontSize', 12);
grid minor
end

if(sel == 0)
figure (3)
plot(XI(1,:)', XI(2,:)', '--k', 'LineWidth',1); 
hold on
plot(xg,yg, 's', 'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6], 'MarkerSize', 10); 
hold on 
plot(x0,y0,'o', 'MarkerEdgeColor','blue','MarkerFaceColor',[.6 .6 1], 'MarkerSize', 10); 
grid minor
title('Trayectoria del robot', 'Interpreter', 'latex', 'Fontsize', 20); 
xlabel('$x$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$y$', 'Interpreter', 'latex', 'Fontsize', 16);
l = legend('Trayectoria','$\zeta_g$','$\zeta_0$','Location', 'best','Orientation', 'vertical');
set(l, 'Interpreter', 'latex', 'FontSize', 12);  
end



% cleanup code goes here: write data to files, etc.
