clc
close all
clear all

%%
% load data
vn310gnss = importgnss("E:\Studium\M2-Inertialnavigation\Uebung\3\Code\vn310-gnss.csv", [2, Inf]);
vn310imu = importimu("E:\Studium\M2-Inertialnavigation\Uebung\3\Code\vn310-imu.csv", [2, Inf]);

% IMU
% acc
accX = vn310imu{:,31};
accY = vn310imu{:,32};
accZ = vn310imu{:,33}; % m/s^2

% dreh
angulX = vn310imu{:,34};
angulY = vn310imu{:,35};
angulZ = vn310imu{:,36}; % rad/s

% time
t_imu = vn310imu{:,3};

% GNSS
% time
t_gnss = vn310gnss{:,3};

% time for kalman filter
t_kalman = zeros(length(t_gnss),1);
for i = 1:length(t_gnss)
    [~,t_kalman(i)] = min(abs(t_gnss(i) - t_imu));
end

% 1
gnss_x1 = vn310gnss{:,31};
gnss_y1 = vn310gnss{:,32};
gnss_z1 = vn310gnss{:,33};

gnss_vx1 = vn310gnss{:,37}; 
gnss_vy1 = vn310gnss{:,38}; 
gnss_vz1 = vn310gnss{:,39}; 

% 2
gnss_x2 = vn310gnss{:,113};

gnss_y2 = vn310gnss{:,114};
gnss_z2 = vn310gnss{:,115};

gnss_vx2 = vn310gnss{:,119}; 
gnss_vy2 = vn310gnss{:,120}; 
gnss_vz2 = vn310gnss{:,121}; 

% richtig internal
gnss_x_inter = vn310gnss{:,76};
gnss_y_inter = vn310gnss{:,77};
gnss_z_inter = vn310gnss{:,78};

% Accuracy
sigma_v = 0.02; % m/s
sigma_xy = 1; % m
sigma_z = 1.5; % m

sigma_acc = 0;
sigma_angul = 0.015; % degree
R = diag([sigma_acc^2,sigma_acc^2,sigma_acc^2,sigma_angul^2,...
    sigma_angul^2,sigma_angul^2]);
%%
% Konstant
f = 1 / 298.256;
e = (2 - f) * f;
a = 6378137;
g_n = [0,0,9.8465]';
tau = 10;           % korrelations Laenge
beta = 1/tau;

w_e = 2 * pi / 86400; % rad/s
omega_ie = [0;0;w_e];
Omega_ie = [0,-w_e,0; w_e,0,0; 0,0,0];

dt = 1/100; % s

acc = [accX,accY,accZ]';
dv = acc * dt;
delta_alpha = [angulX,angulY,angulZ]' * dt; % 

% Transformation
rot_z = [0,-1,0;1,0,0;0,0,1];
rot_x = [1,0,0;0,0,1;0,-1,0];
acc = rot_x * rot_z * acc;
delta_alpha =  rot_x * rot_z * delta_alpha;

% Anfangswert
x_01 = mean(gnss_x1(1:50));
y_01 = mean(gnss_y1(1:50));
z_01 = mean(gnss_z1(1:50));

x_02 = mean(gnss_x2(1:50));
y_02 = mean(gnss_y2(1:50));
z_02 = mean(gnss_z2(1:50));

code = 3;
if code == 1
    x0 = x_01;
    y0 = y_01;
    z0 = z_01;
    hebel_s = [0.0702 0.088 0.474/2]';
elseif code == 2
    x0 = x_02;
    y0 = y_02;
    z0 = z_02;
    hebel_s = [0.0702 0.088 -0.474/2]';
else
    x0 = (x_01+x_02)/2;
    y0 = (y_01+y_02)/2;
    z0 = (z_01+z_02)/2;
    hebel_s = [0.0702 0.088 0]';
end

hebel_b = rot_x * rot_z * hebel_s;

[lon0,lat0,~] = cart2ellip(x0,y0,z0,a,e); % in degree

% degree
Yaw0 = 172;
Pitch0 = 0;
Roll0 = 0;

% initial
len = length(acc);
C_dach = cell(len+1,1);
delta_beta = zeros(3,len);
omega = zeros(3,len);
q = zeros(4,len);
v = zeros(3,len);
pos = zeros(3,len);
pos(:,1) = [x0;y0;z0];
pos(:,2) = [x0;y0;z0];

% now try gnss 1
j = 1;
euler = zeros(3,1);
for i = 2:len
    if i == 2
        C_0 = C_n2e(lat0,lon0) * C_b2n(Yaw0,Pitch0,Roll0);
        C_dach{i-1} = C_0;
        q0 = sqrt(C_0(1,1) + C_0(2,2) + C_0(3,3) + 1) / 2;
        q1 = (C_0(2,3) - C_0(3,2)) / (4 * q0);
        q2 = (C_0(3,1) - C_0(1,3)) / (4 * q0);
        q3 = (C_0(1,2) - C_0(2,1)) / (4 * q0);
        q_0 = [q0 q1 q2 q3]';
        q(:,1) = q_0;
        delta_beta(:,i-1) = delta_alpha(:,i-1) - ...
            inv(C_0) * omega_ie * dt;             % e nach b
        delta_beta(:,i) = delta_alpha(:,i) - ...
            inv(C_dach{i-1}) * omega_ie * dt; 
        omega0 = (3 * delta_beta(:,i-1) - delta_beta(:,i)) / (2 * dt);
        omega(:,i-1) = (delta_beta(:,i-1) + delta_beta(:,i)) / (2 * dt);
        omega(:,i) = (3 * delta_beta(:,i) - delta_beta(:,i-1)) / (2 * dt);
        q(:,i) = RK3(q_0,omega0, omega(:,i-1), omega(:,i),dt);
        q(:,i) = q(:,i) / norm(q(:,i));
        C_dach{i} = DCM(q(:,i));
        v0 = [0;0;0];
        x_0 = [x0;y0;z0];
        g_e =  C_n2e(lat0,lon0) * g_n;
        v(:,i) = v0 + (C_0 * (3 * dv(:,i-1) - dv(:,i)) + ...
                4 * C_dach{i-1} * (dv(:,i-1) + dv(:,i)) + ...
                C_dach{i} * (3 * dv(:,i) - dv(:,i-1))) / 6 ...
                 - (2 * Omega_ie * v0 + Omega_ie * ...
                 Omega_ie * x_0 - g_e) * (2 * dt);
        pos(:,i+1) = pos(:,i) + v(:,i) * dt;
    else
        delta_beta(:,i) = delta_alpha(:,i) - ...
            inv(C_dach{i-1}) * omega_ie * dt; 
        omega(:,i) = (3 * delta_beta(:,i) - delta_beta(:,i-1)) / (2 * dt);
        q(:,i) = RK3(q(:,i-2),omega(:,i-2), omega(:,i-1), omega(:,i),dt);
        q(:,i) = q(:,i) / norm(q(:,i));
        C_dach{i} = DCM(q(:,i));

        [lon,lat,~] = cart2ellip(pos(1,i),pos(2,i),pos(3,i),a,e);
        v(:,i) = v(:,i-2) + (C_dach{i-2} * (3 * dv(:,i-1) - dv(:,i)) + ...
            4 * C_dach{i-1} * (dv(:,i-1) + dv(:,i)) + C_dach{i} * ...
            (3 * dv(:,i) - dv(:,i-1))) / 6 - (2 * Omega_ie * v(:,i-2) ...
            + Omega_ie * Omega_ie * pos(:,i-2) - g_e) * (2 * dt);
        if ismember(i, t_kalman)
            %
            state_var = [ones(3,1)*100^2
                        ones(3,1)*10^2
                        ones(3,1)*0.1^2
                        ones(3,1)*1^2
                        ones(3,1)*0.1^2];
            state_predict = zeros(15,1);
            state_predict(1:3) = state_predict(1:3);
            P = diag(state_var);
            %
            d_pos =  [gnss_x1(j);gnss_y1(j);gnss_z1(j)]...
                - pos(:,i) - hebel_b;
            d_vel = [gnss_vx1(j);gnss_vy1(j);gnss_vz1(j)] - v(:,i);
            j = j+1;
            observation = [d_pos;d_vel];
            % F Matrix
            F = zeros(15);
            F(1:3,4:6) = eye(3);
            F(4:6,1:3) = -Omega_ie * Omega_ie;
            F(4:6,4:6) = -2 * Omega_ie;
            F(4:6,7:9) = [0,-acc(3,i),acc(2,i); 
                          acc(3,i),0,-acc(1,i); 
                          -acc(2,i),acc(1,i),0];
            F(4:6,10:12) = C_dach{i};
            F(7:9,7:9) = -Omega_ie;
            F(7:9,13:15) = -C_dach{i};
            F(10:12,10:12) = -eye(3) * beta;
            F(13:15,13:15) = -eye(3) * beta;
            % G
            G = zeros(15,6);
            G(10:12,1:3) = eye(3) * sqrt(beta * sigma_acc^2);
            G(13:15,4:6) = eye(3) * sqrt(beta * sigma_angul^2);
            %
            H = [eye(6) zeros(6,9)];
            % 
            A = [-F,G * G';zeros(15),F'] * dt;
            B = expm(A);
            B(16:30,16:30);
            B = B';
            Phi = B(16:30,16:30)';
            Q = Phi * B(1:15,16:30);
            % 
            xnnp = Phi * state_predict;
            Pnnp = Phi * P * Phi' + Q;
            K = Pnnp * H' * inv(H * Pnnp * H' + R);
            %
            state_corrected = xnnp + K * (observation - H * xnnp);
            P_1{j} = (eye(15) - K * H) * Pnnp;
            % 
            pos(:,i) = pos(:,i) + state_corrected(1:3);
            v(:,i) = v(:,i) + state_corrected(4:6);
            state_corrected_1(:,j) = state_corrected;
            %
            euler(1) = atan2(C_dach{i}(2,3), ...
                C_dach{i}(3,3))/pi*180; % alpha
            euler(2) = asin(-C_dach{i}(1,3))/pi*180;  % beta
            euler(3) = atan2(C_dach{i}(1,2), ...
                C_dach{i}(1,1))/pi*180; % gamma
            euler = euler + state_corrected(7:9);
            C_dach{i} = rotation(euler(3),'z') * ...
                rotation(euler(2),'x') * rotation(euler(1),'z');
        end
        g_e = C_n2e(lat,lon) * g_n;
        pos(:,i+1) = pos(:,i) + v(:,i) * dt;
    end
end
pos_mitarm = pos;
figure
hold on
plot3(pos(1,:),pos(2,:),pos(3,:))
plot3(gnss_x_inter(201:end),gnss_y_inter(201:end),gnss_z_inter(201:end))
set(gca,'FontSize', 20);
title('mit hebelarm')

%%
% Konstant
f = 1 / 298.256;
e = (2 - f) * f;
a = 6378137;
g_n = [0,0,9.8465]';
tau = 10;           % korrelations Laenge
beta = 1/tau;

w_e = 2 * pi / 86400; % rad/s
omega_ie = [0;0;w_e];
Omega_ie = [0,-w_e,0; w_e,0,0; 0,0,0];

dt = 1/100; % s

acc = [accX,accY,accZ]';
dv = acc * dt;
delta_alpha = [angulX,angulY,angulZ]' * dt; % 

% Transformation
rot_z = [0,-1,0;1,0,0;0,0,1];
rot_x = [1,0,0;0,0,1;0,-1,0];
acc = rot_x * rot_z * acc;
delta_alpha =  rot_x * rot_z * delta_alpha;
% Anfangswert
x_01 = mean(gnss_x1(1:50));
y_01 = mean(gnss_y1(1:50));
z_01 = mean(gnss_z1(1:50));

x_02 = mean(gnss_x2(1:50));
y_02 = mean(gnss_y2(1:50));
z_02 = mean(gnss_z2(1:50));

[lon0,lat0,h0] = cart2ellip(x0,y0,z0,a,e); % in degree

vel_x0 = 0;
vel_y0 = 0;
vel_z0 = 0;

% degree
Yaw0 = 172;
Pitch0 = 0;
Roll0 = 0;

% initial
len = length(acc);
C_dach = cell(len+1,1);
delta_beta = zeros(3,len);
omega = zeros(3,len);
q = zeros(4,len);
v = zeros(3,len);
pos = zeros(3,len);
pos(:,1) = [x0;y0;z0];
pos(:,2) = [x0;y0;z0];

% now try gnss 1
j = 1;
euler = zeros(3,1);
for i = 2:len
    if i == 2
        C_0 = C_n2e(lat0,lon0) * C_b2n(Yaw0,Pitch0,Roll0);
        C_dach{i-1} = C_0;
        q0 = sqrt(C_0(1,1) + C_0(2,2) + C_0(3,3) + 1) / 2;
        q1 = (C_0(2,3) - C_0(3,2)) / (4 * q0);
        q2 = (C_0(3,1) - C_0(1,3)) / (4 * q0);
        q3 = (C_0(1,2) - C_0(2,1)) / (4 * q0);
        q_0 = [q0 q1 q2 q3]';
        q(:,1) = q_0;
        delta_beta(:,i-1) = delta_alpha(:,i-1) - inv(C_0) * omega_ie * dt;             % e nach b
        delta_beta(:,i) = delta_alpha(:,i) - inv(C_dach{i-1}) * omega_ie * dt; 
        omega0 = (3 * delta_beta(:,i-1) - delta_beta(:,i)) / (2 * dt);
        omega(:,i-1) = (delta_beta(:,i-1) + delta_beta(:,i)) / (2 * dt);
        omega(:,i) = (3 * delta_beta(:,i) - delta_beta(:,i-1)) / (2 * dt);
        q(:,i) = RK3(q_0,omega0, omega(:,i-1), omega(:,i),dt);
        q(:,i) = q(:,i) / norm(q(:,i));
        C_dach{i} = DCM(q(:,i));
        v0 = [0;0;0];
        x_0 = [x0;y0;z0];
        g_e =  C_n2e(lat0,lon0) * g_n;
        v(:,i) = v0 + (C_0 * (3 * dv(:,i-1) - dv(:,i)) + ...
                4 * C_dach{i-1} * (dv(:,i-1) + dv(:,i)) + ...
                C_dach{i} * (3 * dv(:,i) - dv(:,i-1))) / 6 ...
                 - (2 * Omega_ie * v0 + Omega_ie * ...
                 Omega_ie * x_0 - g_e) * (2 * dt);
        pos(:,i+1) = pos(:,i) + v(:,i) * dt;
    else
        delta_beta(:,i) = delta_alpha(:,i) - inv(C_dach{i-1})...
            * omega_ie * dt; 
        omega(:,i) = (3 * delta_beta(:,i) - delta_beta(:,i-1))...
            / (2 * dt);
        q(:,i) = RK3(q(:,i-2),omega(:,i-2), omega(:,i-1), omega(:,i),dt);
        q(:,i) = q(:,i) / norm(q(:,i));
        C_dach{i} = DCM(q(:,i));

        [lon,lat,~] = cart2ellip(pos(1,i),pos(2,i),pos(3,i),a,e);
        v(:,i) = v(:,i-2) + (C_dach{i-2} * (3 * dv(:,i-1) - dv(:,i))...
        + 4 * C_dach{i-1} * (dv(:,i-1) + dv(:,i)) + C_dach{i} * ...
        (3 * dv(:,i) - dv(:,i-1))) / 6  - (2 * Omega_ie * v(:,i-2)...
        + Omega_ie * Omega_ie * pos(:,i-2) - g_e) * (2 * dt);
        if ismember(i, t_kalman)
            %
            state_var = [ones(3,1)*100^2
                        ones(3,1)*10^2
                        ones(3,1)*0.1^2
                        ones(3,1)*1^2
                        ones(3,1)*0.1^2];
            state_predict = zeros(15,1);
            P = diag(state_var);
            %
            d_pos =  [gnss_x1(j);gnss_y1(j);gnss_z1(j)] - pos(:,i);
            d_vel = [gnss_vx1(j);gnss_vy1(j);gnss_vz1(j)] - v(:,i);
            j = j+1;
            observation = [d_pos;d_vel];
            % F Matrix
            F = zeros(15);
            F(1:3,4:6) = eye(3);
            F(4:6,1:3) = -Omega_ie * Omega_ie;
            F(4:6,4:6) = -2 * Omega_ie;
            F(4:6,7:9) = [0,-acc(3,i),acc(2,i); 
                          acc(3,i),0,-acc(1,i); 
                          -acc(2,i),acc(1,i),0];
            F(4:6,10:12) = C_dach{i};
            F(7:9,7:9) = -Omega_ie;
            F(7:9,13:15) = -C_dach{i};
            F(10:12,10:12) = -eye(3) * beta;
            F(13:15,13:15) = -eye(3) * beta;
            % G
            G = zeros(15,6);
            G(10:12,1:3) = eye(3) * sqrt(beta * sigma_acc^2);
            G(13:15,4:6) = eye(3) * sqrt(beta * sigma_angul^2);
            %
            H = [eye(6) zeros(6,9)];
            % 
            A = [-F,G * G';zeros(15),F'] * dt;
            B = expm(A);
            B(16:30,16:30);
            B = B';
            Phi = B(16:30,16:30)';
            Q = Phi * B(1:15,16:30);
            % 
            xnnp = Phi * state_predict;
            Pnnp = Phi * P * Phi' + Q;
            K = Pnnp * H' * inv(H * Pnnp * H' + R);
            %
            state_corrected = xnnp + K * (observation - H * xnnp);
            P_2{j} = (eye(15) - K * H) * Pnnp;
            % 
            pos(:,i) = pos(:,i) + state_corrected(1:3);
            v(:,i) = v(:,i) + state_corrected(4:6);
            state_corrected_2(:,j) = state_corrected;
            %
            euler(1) = atan2(C_dach{i}(2,3), ...
                C_dach{i}(3,3))/pi*180; % alpha
            euler(2) = asin(-C_dach{i}(1,3))/pi*180;  % beta
            euler(3) = atan2(C_dach{i}(1,2), ...
                C_dach{i}(1,1))/pi*180; % gamma
            euler = euler + state_corrected(7:9);
            C_dach{i} = rotation(euler(3),'z') * ...
                rotation(euler(2),'x') * rotation(euler(1),'z');
        end
        g_e = C_n2e(lat,lon) * g_n;
        pos(:,i+1) = pos(:,i) + v(:,i) * dt;
    end
end
pos_ohnearm = pos;
figure
hold on
plot3(pos(1,:),pos(2,:),pos(3,:))
plot3(gnss_x_inter(201:end),gnss_y_inter(201:end),gnss_z_inter(201:end))
set(gca,'FontSize', 20);
title('ohne hebelarm')

diff = pos_ohnearm - pos_mitarm;
diff2 = state_corrected_1 - state_corrected_2;
figure
subplot(3,1,1)
plot(pos_ohnearm(1,100:end) - pos_mitarm(1,100:end))
title('x')
set(gca,'FontSize', 20);
subplot(3,1,2)
plot(pos_ohnearm(2,100:end) - pos_mitarm(2,100:end))
title('y')
set(gca,'FontSize', 20);
subplot(3,1,3)
plot(pos_ohnearm(3,100:end) - pos_mitarm(3,100:end))
title('z')
set(gca,'FontSize', 20);

pos_compare_mit = pos_mitarm(:,t_kalman); 
pos_compare_ohne = pos_ohnearm(:,t_kalman); 

f_mit = pos_compare_mit - [gnss_x_inter,gnss_y_inter,gnss_z_inter]';
f_ohne = pos_compare_ohne - [gnss_x_inter,gnss_y_inter,gnss_z_inter]';

figure
subplot(3,1,1)
plot(f_mit(1,201:end))
title('x')
set(gca,'FontSize', 20);
subplot(3,1,2)
plot(f_mit(2,201:end))
title('y')
set(gca,'FontSize', 20);
subplot(3,1,3)
plot(f_mit(3,201:end))
title('z')
set(gca,'FontSize', 20);

% figure
% subplot(3,1,1)
% plot(f_ohne(1,1:200)-mean(f_ohne(1,1:200)))
% title('x')
% set(gca,'FontSize', 20);
% subplot(3,1,2)
% plot(f_ohne(2,1:200)-mean(f_ohne(2,1:200)))
% title('y')
% set(gca,'FontSize', 20);
% subplot(3,1,3)
% plot(f_ohne(3,1:200)-mean(f_ohne(3,1:200)))
% title('z')
% set(gca,'FontSize', 20);
