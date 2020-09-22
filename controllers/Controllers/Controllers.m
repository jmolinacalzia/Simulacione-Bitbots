% MATLAB controller for Webots
% File:          Controllers.m
% Date:
% Description:
% Author:
% Modifications:

% uncomment the next two lines if you want to use
% MATLAB's desktop to interact with the controller:
desktop;
%keyboard;



TIME_STEP = 50; 
left_wheel = wb_robot_get_device('left wheel');
right_wheel = wb_robot_get_device('right wheel');
position_sensor = wb_robot_get_device('gps');
orientation_sensor = wb_robot_get_device('compass');



wb_gps_enable(position_sensor, TIME_STEP);
wb_compass_enable(orientation_sensor, TIME_STEP);


wb_motor_set_position(left_wheel, inf);
wb_motor_set_position(right_wheel, inf);
wb_motor_set_velocity(left_wheel, 0.0);
wb_motor_set_velocity(right_wheel, 0.0);
wb_gps_enable(position_sensor, TIME_STEP);
wb_compass_enable(orientation_sensor, TIME_STEP);


MAX_SPEED = 12.3; 
WHEEL_RADIUS = 195/2000; 
DISTANCE_CENTER = 381/2000;
ell = 0.1; 
pos0 = wb_gps_get_values(position_sensor); 
norte = wb_compass_get_values(orientation_sensor);
theta0 = atan2(norte(1), norte(3));

%================================================================

% Variables de control

sel = 'fg'; %--> fg, traj
ctrl = 'lqi'; % --> pid, pidexp, nonlin, lqr, lqi

%================================================================
  
% Metas: 

if (strcmp(sel,'fg'))

  xg = 7; 
  yg = 5;
  thetag = pi; 
  
elseif (strcmp(sel,'traj'))
  yg = 0; %definir trayectoria aqui
end 
%===============================================================
% Difeomorfismo para LQR y LQI

difeo = @(theta,mu) [1,0;0,1/ell]*[cos(theta), -sin(theta); sin(theta), cos(theta)]'*[mu(1); mu(2)]; 


%================================================================

% Condiciones iniciales 

pos0 = wb_gps_get_values(position_sensor); 
pos = pos0'; 
norte = wb_compass_get_values(orientation_sensor);
theta0 = atan2(norte(1), norte(3));
theta = theta0; 

xi = [pos(1);pos(3);theta]; 

u0 = [0;0]; 

%XI = xi; 
%U = u0; 


dt = TIME_STEP/1000; 

%==================================================================

% Controladores

%% PID 
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

%% Controlador no lineal 
k_rho = 10; 
k_alpha = 20; 
k_beta = -10; 

%% LQR
Klqr = [1,0;0,1]; 
mu = [0;0]

%% LQI
Klqi = [1,0,3.1623,0;0,1,0,3.1623]; 
sigma = [0;0]; 
mu = [0;0]; 
Ekx = 0; ekx = 0; 
Eky = 0; eky = 0; 
v = 0; 
w = 0; 


%==================================================================

% Parámetros de simulación 
t0 = 0; 
tf = 40; 
count = 1; 
tmax = (tf-t0)/(TIME_STEP/1000); 


% Históricos para gráficas 
U = u0; 
XI = xi; 
PHIL = 0; 
PHIR = 0; 


%==================================================================
 
% get and enable devices, e.g.:
%  camera = wb_robot_get_device('camera');
%  wb_camera_enable(camera, TIME_STEP);
%  motor = wb_robot_get_device('motor');

% main loop:
% perform simulation steps of TIME_STEP milliseconds
% and leave the loop when Webots signals the termination
%
while wb_robot_step(TIME_STEP) ~= -1

pos = wb_gps_get_values(position_sensor); 
pos = pos'; 
norte = wb_compass_get_values(orientation_sensor);
theta = atan2(norte(1), norte(3));
  
  xi = [pos(1);pos(3);theta]; 

  
  if(strcmp(ctrl, 'pid'))
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
       constP = MAX_SPEED * (1-exp(-alpha*eP^2)) / eP;
       v = constP*eP;
       
       % Ori --> Velocidad angular 
       eO_D = eO - eO_1;
       EO = EO + eO;
       w = kpO*eO + kiO*EO + kdO*eO_D;
       eO_1 = eO;       
       
       % Vector de entradas
       u = [v; w];
       
  elseif(strcmp(ctrl, 'lqi'))

       x = xi(1:2); 
       ref = [xg;yg]; 
        
        
       ekx = x(1) - ref(1); 
       Ekx = Ekx+ekx; 
        
       eky = x(2) - ref(2); 
       Eky = Eky+eky; 
        
        
       mu(1) = -5*sqrt(3)*x(1) - 1.5*Ekx*dt; 
       mu(2) = -5*sqrt(3)*x(2) - 1.5*Eky*dt; 
      
      
      
        % Vector de entradas
      
       u = -difeo(theta,mu);
       v = u(1); w = u(2); 
         
        
       wb_console_print(sprintf('v = %2.2f\n w = %2.2f', v,w),WB_STDOUT); 
       wb_console_print(sprintf('theta = %2.2f\n' , theta),WB_STDOUT);
  elseif(strcmp(ctrl, 'lqr'))
  
       x = xi(1); 
       y = xi(2);
       e = [x - xg; y - yg];
         
       mu(1) = -3*e(1); 
       mu(2) = -3*e(2);  

       % Vector de entradas 
       u = -difeo(theta,mu); 
       v = u(1); w = u(2); 
       
       wb_console_print(sprintf('v = %2.2f\n w = %2.2f', v, w),WB_STDOUT); 
  elseif(strcmp(ctrl, 'nonlin'))
       x = xi(1); 
       y = xi(2); 
       theta = xi(3);
       e = [xg - x; yg - y];
            
       rho = norm(e);
       alpha = -theta + atan2(e(2), e(1));
       alpha = atan(sin(alpha)/cos(alpha)); 
       beta = -theta - alpha;
            
       v = k_rho*rho;
       w = k_alpha*alpha + k_beta*beta;
            
       % Vector de entradas
       u = [v;w]; 
  end
  % read the sensors, e.g.:
  %  rgb = wb_camera_get_image(camera);

  % Process here sensor data, images, etc.

  % send actuator commands, e.g.:
  %  wb_motor_set_postion(motor, 10.0);

  % if your code plots some graphics, it needs to flushed like this:
  
  
  drawnow;
  speedL = (v-w*DISTANCE_CENTER/WHEEL_RADIUS);
  speedR = (v+w*DISTANCE_CENTER/WHEEL_RADIUS);   
  
  
  if(strcmp(ctrl,'nonlin'))
  wb_motor_set_velocity(left_wheel, speedR); 
  wb_motor_set_velocity(right_wheel, speedL);
  
  elseif(strcmp(ctrl,'lqi'))
  wb_motor_set_velocity(left_wheel, speedL/50); 
  wb_motor_set_velocity(right_wheel, speedR/50);
  else
  wb_motor_set_velocity(left_wheel, speedL); 
  wb_motor_set_velocity(right_wheel, speedR);

  end
  
 

  wb_console_print(sprintf('SpeedL = %2.2f\n SpeedR = %2.2f', speedL, speedR),WB_STDOUT); 
  
  U = [U,u]; 
  XI = [XI,xi]; 
  PHIL = [PHIL,speedL]; 
  PHIR = [PHIR, speedR]; 
  
  count = count + 1; 
  if (count > tmax)
  break; 
  end
  
end



t = linspace(t0,tf,count); 

figure(1); 
plot(t, XI', 'LineWidth', 1); 
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\mathbf{x}(t)$', 'Interpreter', 'latex', 'Fontsize', 16);
hold on; 
yline(yg, '--', '$y_g$','Interpreter', 'latex'); 
yline(xg, '-.', '$x_g$','Interpreter', 'latex'); 
l = legend('$\zeta(t)$', '$y(t)$', '$\theta(t)$', '$y_g$', '$x_g$',...
 'Location', 'best', ...
    'Orientation', 'vertical');
set(l, 'Interpreter', 'latex', 'FontSize', 12);
title('Salida del controlador', 'Interpreter', 'latex', 'Fontsize', 20); 
ylim([-2.5,8.5])
grid minor;




figure(2); 
plot(t,U', 'LineWidth', 1); 
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\mathbf{u}(t)$', 'Interpreter', 'latex', 'Fontsize', 16);
l = legend('$v(t)$', '$\omega(t)$',...
 'Location', 'best', ...
    'Orientation', 'vertical');
set(l, 'Interpreter', 'latex', 'FontSize', 12);
grid minor; 


figure(3); 
plot(t,PHIL', 'LineWidth', 1); 
hold on; 
plot(t,PHIR', 'LineWidth', 1); 
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', 16);
ylabel('$\mathbf{\phi}(t)$', 'Interpreter', 'latex', 'Fontsize', 16);
l = legend('$\phi_{L}(t)$', '$\phi_{R}(t)$',...
 'Location', 'best', ...
    'Orientation', 'vertical');
set(l, 'Interpreter', 'latex', 'FontSize', 12);
grid minor; 
