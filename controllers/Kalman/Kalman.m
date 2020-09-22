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

x0 = -5 + delta0; 
y0 = -5 + delta0; 
theta0 = pi; 


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

% Parámetros de simulación 
t0 = 0; 
tf = 25; 
count = 1; 
tmax = (tf-t0)/(TIME_STEP/1000); 



cam_count = 0; 
y_1 = 0; 


% Loop de simulación 
while wb_robot_step(TIME_STEP) ~= -1

  % Controlador 
  

  phiR = 10; 
  phiL = 12; 
  
  mu = [phiR; phiL]; 
  
  wb_motor_set_velocity(left_wheel, phiL); 
  wb_motor_set_velocity(right_wheel, phiR);

  
  
  
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
  fprintf("Corrección\n")
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
  % Guardar salida del EKF
  XHAT = [XHAT, xhat_post]; 
  
  % Trayectorias de estado reales 
  XI = [XI,xi]; 
  MU = [MU,mu]; 
  
  % Mediciones de sensores
  
  CAM = [CAM, CAM_SENSOR]; 
  ENCR = [ENCR, encr_pulses]; 
  ENCL = [ENCL, encl_pulses]; 
  
  
    
  
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
l = legend('$x(t)$', '$y(t)$', '$\theta(t)$','$\hat{x}(t)$','$\hat{y}(t)$','$\hat{\theta}(t)$',...
 'Location', 'best', ...
    'Orientation', 'vertical');
set(l, 'Interpreter', 'latex', 'FontSize', 12);
grid minor;


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


% cleanup code goes here: write data to files, etc.
