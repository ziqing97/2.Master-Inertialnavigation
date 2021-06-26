close all;
format longg;
clear;
clc;

% Autor: Carsten Helfert, 3318553

%% INav Uebung 2

% Aufgabe 1 a)
% Aufgabe 3 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Nr3 = true;         % Nr3 = true --> Ergebnisse Aufgabe 3
                    % Nr3 = false --> Ergebnisse Aufgabe 2

update_rate = 5;
% update_rate = 60;
% update_rate = 120;    % Update alle update_rate [s] (siehe Nr3 a)


data = table2array(readtable('vn-data-static.csv'));

% Acc 28-30
acc = data(:,28:30);

% Gyro 31-33
gyro = data(:,31:33);

% Time [ns]
time = data(:,4)*1e-9;  % [s]

% Datenlücken numerisch suchen
acc_idx = find(isnan(acc))';
gyro_idx = find(isnan(gyro))';

% Datenlücken graphisch ermitteln
% Beschleunigung darstellen
figure;
subplot(3,1,1);
hold on;
title('Beschleunigung in der x-Komponente');
plot(time, acc(:,1), '.', 'MarkerSize', 0.1);
ylabel('Beschleunigung $\left[\frac{m}{s^2}\right]$','Interpreter','latex')
xlabel('Zeit [ns]')
hold off;

subplot(3,1,2);
hold on;
title('Beschleunigung in der y-Komponente');
plot(time, acc(:,2), '.', 'MarkerSize', 0.1);
ylabel('Beschleunigung $\left[\frac{m}{s^2}\right]$','Interpreter','latex')
xlabel('Zeit [s]')
hold off;

subplot(3,1,3);
hold on;
plot(time, acc(:,3), '.', 'MarkerSize', 0.1);
title('Beschleunigung in der z-Komponente');
ylabel('Beschleunigung $\left[\frac{m}{s^2}\right]$','Interpreter','latex')
xlabel('Zeit [s]')
hold off;
if ~isfile('acc.png')
saveas(gcf,'acc.png');
end

% Drehraten darstellen
figure;
subplot(3,1,1);
hold on;
title('Drehraten in der x-Komponente');
plot(time, gyro(:,1), '.', 'MarkerSize', 0.1);
ylabel('Drehraten $\left[\frac{rad}{s}\right]$','Interpreter','latex')
xlabel('Zeit [s]')
hold off;

subplot(3,1,2);
hold on;
title('Drehraten in der y-Komponente');
plot(time, gyro(:,2), '.', 'MarkerSize', 0.1);
ylabel('Drehraten $\left[\frac{rad}{s}\right]$','Interpreter','latex')
xlabel('Zeit [s]')
hold off;

subplot(3,1,3);
hold on;
title('Drehraten in der z-Komponente');
plot(time, gyro(:,3), '.', 'MarkerSize', 0.1);
ylabel('Drehraten $\left[\frac{rad}{s}\right]$','Interpreter','latex')
xlabel('Zeit [s]')
hold off;
if ~isfile('gyro.png')
saveas(gcf,'gyro.png');
end

% Berechnen der Standardabweichungen
sigma_acc_x = sqrt(var(acc(:,1)));
sigma_acc_y = sqrt(var(acc(:,2)));
sigma_acc_z = sqrt(var(acc(:,3)));
sigma_acc_pos = sqrt(sigma_acc_x^2+sigma_acc_y^2+sigma_acc_z^2);

sigma_gyro_x = sqrt(var(gyro(:,1)));
sigma_gyro_y = sqrt(var(gyro(:,2)));
sigma_gyro_z = sqrt(var(gyro(:,3)));
sigma_gyro_pos = sqrt(sigma_gyro_x^2+sigma_gyro_y^2+sigma_gyro_z^2);


%% Aufgabe 2

% Winkel omega bestimmen
% Initialwerte für Yaw, Pitch, Roll
delta_t = 0.02; % [s]
d2r = pi/180;
Y0 = -118*d2r;  % [rad]
P0 = 0; % horizontierte IMU
R0 = 0; % -> in lokaler Tangentialebene
p_0 = zeros(6, length(gyro));

% Startwerte der Quaternionen bestimmen
lat = 48.78070192*d2r; % [rad]
lon = 9.17158708*d2r; % [rad]
alt = 326.568; % [m]
pos = geod2cart(lat, lon, alt);
pos = [pos pos zeros(3,length(gyro)-1)];
vel = zeros(3,length(gyro)+1);
euler = zeros(3,length(gyro));
omega_ie = [0 0 2*pi/86400]';  % [rad/s]
g_n = [0 0 9.8465]';    % [m/s^2]
Omega_ie = [0 -2*pi/86400 0
            2*pi/86400 0 0
            0 0 0];

% Initialisierung der DCM_pe
DCM_bn = [cos(Y0)*cos(P0) cos(Y0)*sin(P0)*sin(R0)-sin(Y0)*cos(R0) cos(Y0)*sin(P0)*cos(R0)+sin(Y0)*sin(R0);
         sin(Y0)*cos(P0) sin(Y0)*sin(P0)*sin(R0)+cos(Y0)*cos(R0), sin(Y0)*sin(P0)*cos(R0)-cos(Y0)*sin(R0);
         -sin(P0) cos(P0)*sin(R0)  cos(P0)*cos(R0)];
DCM_ne = [-sin(lat)*cos(lon) -sin(lon) -cos(lat)*cos(lon)
          -sin(lat)*sin(lon) cos(lon) -cos(lat)*sin(lon)
          cos(lat) 0 -sin(lat)];
DCM_0 = DCM_ne * DCM_bn;
DCM_1 = DCM_0;

euler(1,1) = atan2(DCM_0(2,3), DCM_0(3,3))/d2r;
euler(2,1) = asin(-DCM_0(1,3))/d2r;
euler(3,1) = atan2(DCM_0(1,2), DCM_0(1,1))/d2r;

% Extrahieren der Quaternionen aus der DCM
q0 = sqrt(DCM_0(1,1)+DCM_0(2,2)+DCM_0(3,3)+1)/2;
q1 = (DCM_0(2,3)-DCM_0(3,2))/(4*q0);
q2 = (DCM_0(3,1)-DCM_0(1,3))/(4*q0);
q3 = (DCM_0(1,2)-DCM_0(2,1))/(4*q0);
q_0 = [q0 q1 q2 q3]';
q_1 = q_0;

% Aufgabe 3 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (Nr3)
    % Aufgabe 3 - Kalman Filter - Initialwerte
    x_var = [ones(3,1)*100^2
             ones(3,1)*10^2
             ones(3,1)*0.1^2
             ones(3,1)*1^2
             ones(3,1)*0.1^2];
    var_pos = 0.001^2;	% [m^2]
    var_vel = 0.001^2;  % [m^2/s^2]
    var_acc = 1e-5;     % [m^2/s^4]
    var_gyro = 1e-8;    % [rad^2/s^2]
    tau = 10;
    P = {};
    P{1} = eye(15).*x_var';
    x = zeros(15,1);
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
j = 1;
for i=1:length(gyro)-1
    if (i ~= 1)
        DCM_0 = DCM_1;
        DCM_1 = DCM_2;
        q_0 = q_1;
        q_1 = q_2;
        [lat, lon, alt] = cart2geod(pos(:,i-1));
        DCM_ne = [-sin(lat)*cos(lon) -sin(lon) -cos(lat)*cos(lon)
              -sin(lat)*sin(lon) cos(lon) -cos(lat)*sin(lon)
              cos(lat) 0 -sin(lat)];
    end
    
    % Berechnen der Euler-Winkel
    euler(1,i+1) = atan2(DCM_1(2,3), DCM_0(3,3))/d2r;
    euler(2,i+1) = asin(-DCM_1(1,3))/d2r;
    euler(3,i+1) = atan2(DCM_1(1,2), DCM_0(1,1))/d2r;
    
    % Korrigieren der Orientierung
    delta_beta_0 = gyro(i,:)' - DCM_0' * omega_ie * delta_t;
    delta_beta_1 = gyro(i+1,:)' - DCM_1' * omega_ie * delta_t;

    omega_0 = (3*delta_beta_0 - delta_beta_1) / (2*delta_t);
    omega_1 = (delta_beta_0 + delta_beta_1) / (2*delta_t);
    omega_2 = (3*delta_beta_1 - delta_beta_0) / (2*delta_t);

    % Berechnen der Quaternionen mittels Runge-Kutta
    q_2 = RK3(@odefun, q_0, delta_t, omega_0, omega_1, omega_2);

    % Berechnen der neuen DCM
    DCM_2 = [q_2(1)^2+q_2(2)^2-q_2(3)^2-q_2(4)^2 2*(q_2(2)*q_2(3)+q_2(4)*q_2(1)) 2*(q_2(2)*q_2(4)-q_2(3)*q_2(1))
             2*(q_2(2)*q_2(3)-q_2(4)*q_2(1)) q_2(1)^2-q_2(2)^2+q_2(3)^2-q_2(4)^2 2*(q_2(3)*q_2(4)+q_2(2)*q_2(1))
             2*(q_2(2)*q_2(4)+q_2(3)*q_2(1)) 2*(q_2(3)*q_2(4)-q_2(2)*q_2(1)) q_2(1)^2-q_2(2)^2-q_2(3)^2+q_2(4)^2];
    
    % Aufgabe 3 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (Nr3)
        % Übergabeparameter
        x_0 = x;
        pos_error = pos(:,1) - pos(:,i+1);
        vel_error = vel(:,1) - vel(:,i+1);
        acc_val = acc(i+1,:)';
        Cpe = DCM_2;
        
        % Kovarianzmatrix

        
        % Kalman Filterung
        [x, P{j+1}] = Kalman_Filter(x_0, pos_error, var_pos, vel_error, ...
                               var_vel, P{j}, var_acc, var_gyro, ...
                               acc_val, Cpe, Omega_ie, tau, delta_t);
       
        if(round(mod(i-2,update_rate/delta_t)) == 0) ...
        %|| round(mod(i-2,update_rate/delta_t)) == 1)
            p_0(1,j) = P{j}(1,1);
            p_0(2,j) = P{j}(2,2);
            p_0(3,j) = P{j}(3,3);
            p_0(4,j) = P{j}(4,4);
            p_0(5,j) = P{j}(5,5);
            p_0(6,j) = P{j}(6,6);
            [x, P{j+1}] = Kalman_Filter(x_0, pos_error, var_pos, vel_error, ...
                               var_vel, P{j}, var_acc, var_gyro, ...
                               acc_val, Cpe, Omega_ie, tau, delta_t);
           j = j+1;
            % Korrektur mit x alle "update_rate" [s]
            % Update für P jede Epoche
            pos(:,i+1) = pos(:,i+1) + x(1:3);
            vel(:,i+1) = vel(:,i+1) + x(4:6);
            x_0 = [zeros(6,1)
                   x_0(7:end)];
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         
    % Integration der Positions- und Geschwindigkeitsgleichungen
    g_e = DCM_ne * g_n;
    vel(:,i+2) = vel(:,i) ...
                 + (DCM_0*(3*acc(i,:)'*delta_t-acc(i+1,:)'*delta_t)+4*DCM_1*(acc(i,:)'*delta_t+acc(i+1,:)'*delta_t)+...
                    DCM_2*(3*acc(i+1,:)'*delta_t-acc(i,:)'*delta_t))/6 ...
                 - (2*Omega_ie*vel(:,i)+Omega_ie*Omega_ie*pos(:,i)-g_e)*2*delta_t;
    pos(:,i+2) = pos(:,i+1) ...
                 + vel(:,i+1)*delta_t;
end

if (Nr3 == false)
    figure;
    hold on;
    title('Berechnete Positionen');
    plot3(pos(1,:), pos(2,:), pos(3,:));
    ylabel('X-Koordinate [m]');
    xlabel('Y-Koordinate [m]');
    zlabel('Z-Koordinate [m]');
    hold off;
    if ~isfile('position.png')
    saveas(gcf,'position.png');
    end

    figure;
    subplot(3,1,1);
    hold on;
    title('Berechnete Geschwindigkeiten');
    plot(time, vel(1,2:end));
    hold off;

    subplot(3,1,2);
    hold on;
    ylabel('Geschwindigkeit $\left[\frac{m}{s}\right]$','Interpreter','latex','FontSize',12)
    plot(time, vel(2,2:end));
    hold off;

    subplot(3,1,3);
    hold on;
    plot(time, vel(3,2:end));

    xlabel('Zeit [s]','FontSize',12)
    hold off;
    if ~isfile('velocity.png')
    saveas(gcf,'velocity.png');
    end

    figure;
    subplot(3,1,1);
    hold on;
    title('Euler-Winkel');
    plot(time, euler(1,:));
    ylabel('\alpha [°]','FontSize',12);
    hold off;

    subplot(3,1,2);
    hold on;
    plot(time, euler(2,:));
    ylabel('\beta [°]','FontSize',12);
    hold off;

    subplot(3,1,3);
    hold on;
    plot(time, euler(3,:));
    ylabel('\gamma [°]','FontSize',12);
    xlabel('Zeit [s]','FontSize',12);
    hold off;
    if ~isfile('euler.png')
    saveas(gcf,'euler.png');
    end
end

% Aufgabe 3 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (Nr3)
    figure;
    hold on;
    title('Berechnete Positionen');
    plot3(pos(1,:), pos(2,:), pos(3,:));
    ylabel('X-Koordinate [m]');
    xlabel('Y-Koordinate [m]');
    zlabel('Z-Koordinate [m]');
    hold off;
    if ~isfile('position_nr3.png')
    saveas(gcf,'position_nr3.png');
    end
    
    figure;
    hold on;
    title('Berechnete Positionen mit reduzierten Koordinaten');
    plot3(pos(1,:) - pos(1,1), pos(2,:) - pos(2,1), pos(3,:) - pos(3,1));
    ylabel('Reduzierte X-Koordinate [m]');
    xlabel('Reduzierte Y-Koordinate [m]');
    zlabel('Reduzierte Z-Koordinate [m]');
    hold off;
    if ~isfile('position_nr3_red.png')
    saveas(gcf,'position_nr3_red.png');
    end

    figure;
    subplot(3,1,1);
    hold on;
    title('Berechnete Geschwindigkeiten');
    plot(time, vel(1,2:end));
    hold off;

    subplot(3,1,2);
    hold on;
    ylabel('Geschwindigkeit $\left[\frac{m}{s}\right]$','Interpreter','latex','FontSize',12)
    plot(time, vel(2,2:end));
    hold off;

    subplot(3,1,3);
    hold on;
    plot(time, vel(3,2:end));

    xlabel('Zeit [s]','FontSize',12)
    hold off;
    if ~isfile('velocity_nr3.png')
    saveas(gcf,'velocity_nr3.png');
    end
    
    figure;
    hold on;
    title('Betrag der berechneten Geschwindigkeiten');
    plot(time, sqrt(vel(1,2:end).^2+vel(2,2:end).^2+vel(3,2:end).^2));
    hold off;
    ylabel('Geschwindigkeit $\left[\frac{m}{s}\right]$','Interpreter','latex')
    xlabel('Zeit [s]')
    hold off;
    if ~isfile('velocity_nr3_abs.png')
    saveas(gcf,'velocity_nr3_abs.png');
    end
    
    p0id = find(p_0(1,:) ~= 0);
    p_0 = p_0(:,p0id);
    figure;
    subplot(3,1,1);
    hold on;
    title('Varianzen der Geschwindigkeit');
    plot(p_0(1,3:end));
    hold off;
    ylabel('Varianx der x-Position [m]')
    hold off;
    subplot(3,1,2);
    hold on;
    title('Varianzen der Geschwindigkeit');
    plot(p_0(2,3:end));
    hold off;
    ylabel('Varianx der y-Position [m]')
    hold off;
    subplot(3,1,3);
    hold on;
    title('Varianzen der Geschwindigkeit');
    plot(p_0(3,3:end));
    hold off;
    ylabel('Varianx der z-Position [m]')
    hold off;
    if ~isfile('var_pos.png')
    saveas(gcf,'var_pos.png');
    end
    
    figure;
    subplot(3,1,1);
    hold on;
    title('Varianzen der Position');
    plot(p_0(4,3:end));
    hold off;
    ylabel('Varianx der x-Geschwindigkeit [m]')
    hold off;
    subplot(3,1,2);
    hold on;
    title('Varianzen der Position');
    plot(p_0(5,3:end));
    hold off;
    ylabel('Varianx der y-Geschwindigkeit [m]')
    hold off;
    subplot(3,1,3);
    hold on;
    title('Varianzen der Position');
    plot(p_0(6,3:end));
    hold off;
    ylabel('Varianx der z-Geschwindigkeit [m]')
    hold off;
    if ~isfile('variance_vel.png')
    saveas(gcf,'variance_vel.png');
    end
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

figure
subplot(3,1,1)
plot(time,pos(1,2:end) - pos(1,1));
xlabel('Zeit (s)')
ylabel('x Abweichung (m)')
subplot(3,1,2)
plot(time,pos(2,2:end) - pos(2,1));
xlabel('Zeit (s)')
ylabel('y Abweichung (m)')
subplot(3,1,3)
plot(time,pos(3,2:end) - pos(3,1));
xlabel('Zeit (s)')
ylabel('z Abweichung (m)')

