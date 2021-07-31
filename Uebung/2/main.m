%% INav Uebung 2
% Initial
clc
close all
clearvars

%% data
vndatastatic = importfile("E:\Studium\M2-Inertialnavigation\Uebung\2\vn-data-static.csv", [2, Inf]);
%
w_e = 2 * pi / 86400; % rad/s
g_n = [0,0,9.8465]';

t = vndatastatic{:,4} * 1e-9; % s

accX = vndatastatic{:,28};
accY = vndatastatic{:,29};
accZ = vndatastatic{:,30}; % m/s^2

std_accX = std(accX);
std_accY = std(accY);
std_accZ = std(accZ);

figure
subplot(3,1,1)
plot(t,accX)
xlabel("time (ns)")
ylabel("AccX (rad/s)")
subplot(3,1,2)
plot(t,accY)
xlabel("time (ns)")
ylabel("AccY (rad/s)")
subplot(3,1,3)
plot(t,accZ)
xlabel("time (ns)")
ylabel("AccZ (rad/s)")


GyroX = vndatastatic{:,31};
GyroY = vndatastatic{:,32};
GyroZ = vndatastatic{:,33}; % rad/s

std_GyroX = std(accX);
std_GyroY = std(accY);
std_GyroZ = std(accZ);

figure
subplot(3,1,1)
plot(t,GyroX)
xlabel("time (ns)")
ylabel("GyroX (rad/s)")
subplot(3,1,2)
plot(t,GyroY)
xlabel("time (ns)")
ylabel("GyroY (rad/s)")
subplot(3,1,3)
plot(t,GyroZ)
xlabel("time (ns)")
ylabel("GyroZ (rad/s)")

%%
len = length(accX);
dt = 1 / 50;

% Anfangswert
lat0 = 48.78070192;
lon0 = 9.17158708;
h0 = 326.568;
lat0 = lat0 / 180 * pi;
lon0 = lon0 / 180 * pi;

pitch0 = 0;
roll0 = 0;
yaw0 = -118;
yaw0 = yaw0 /180 * pi;

f = 1 / 298.256;
e = (2 - f) * f;
a = 6378137;
[x0, y0, z0]  = ellip2cart(lon0,lat0,h0,a,e);

omega_ie = [0;0;w_e];
Omega_ie = [0,-w_e,0; w_e,0,0; 0,0,0];

delta_alpha = [GyroX,GyroY,GyroZ]';
acc = [accX,accY,accZ]';
dv = acc * dt;

C = cell(len+1,1);
C_dach = cell(len+1,1);
v = zeros(3,len);
delta_beta = zeros(3,len);
omega = zeros(3,len);
pos = zeros(3,len);
pos(:,1) = [x0;y0;z0];
pos(:,2) = [x0;y0;z0];
q = zeros(4,len);
euler = zeros(3,len-1);

for i = 2:len
    if i == 2
        C_0 = C_n2e(lat0,lon0) * C_b2n(yaw0,pitch0,roll0);
        C_dach{i-1} = C_0;
        euler(1,i-1) = atan2(C_dach{i-1}(2,3), C_dach{i-1}(3,3))/pi*180;
        euler(2,i-1) = asin(-C_dach{i-1}(1,3))/pi*180;
        euler(3,i-1) = atan2(C_dach{i-1}(1,2), C_dach{i-1}(1,1))/pi*180;
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
        C_dach{i} = [q(1,i)^2+q(2,i)^2-q(3,i)^2-q(4,i)^2 2*(q(2,i)*q(3,i)+q(4,i)*q(1,i)) 2*(q(2,i)*q(4,i)-q(3,i)*q(1,i))
             2*(q(2,i)*q(3,i)-q(4,i)*q(1,i)) q(1,i)^2-q(2,i)^2+q(3,i)^2-q(4,i)^2 2*(q(3,i)*q(4,i)+q(2,i)*q(1,i))
             2*(q(2,i)*q(4,i)+q(3,i)*q(1,i)) 2*(q(3,i)*q(4,i)-q(2,i)*q(1,i)) q(1,i)^2-q(2,i)^2-q(3,i)^2+q(4,i)^2];
        v0 = [0;0;0];
        x_0 = [x0;y0;z0];
        g_e =  C_n2e(lat0,lon0) * g_n;
        v(:,i) = v0 + (C_0 * (3 * dv(:,i-1) - dv(:,i)) + ...
                4 * C_dach{i-1} * (dv(:,i-1) + dv(:,i)) + ...
                C_dach{i} * (3 * dv(:,i) - dv(:,i-1))) / 6 ...
                 - (2 * Omega_ie * v0 + Omega_ie * Omega_ie * x_0 - g_e) * (2 * dt);
        pos(:,i+1) = pos(:,i) + v(:,i) * dt;
    else
        delta_beta(:,i) = delta_alpha(:,i) - inv(C_dach{i-1}) * omega_ie * dt; 
        omega(:,i) = (3 * delta_beta(:,i) - delta_beta(:,i-1)) / (2 * dt);
        q(:,i) = RK3(q(:,i-2),omega(:,i-2), omega(:,i-1), omega(:,i),dt);
        q(:,i) = RK3_2(@odefun, q(:,i-2), dt, omega(:,i-2), omega(:,i-1), omega(:,i));
        q(:,i) = q(:,i) / norm(q(:,i));
        C_dach{i} = [q(1,i)^2+q(2,i)^2-q(3,i)^2-q(4,i)^2 2*(q(2,i)*q(3,i)+q(4,i)*q(1,i)) 2*(q(2,i)*q(4,i)-q(3,i)*q(1,i))
             2*(q(2,i)*q(3,i)-q(4,i)*q(1,i)) q(1,i)^2-q(2,i)^2+q(3,i)^2-q(4,i)^2 2*(q(3,i)*q(4,i)+q(2,i)*q(1,i))
             2*(q(2,i)*q(4,i)+q(3,i)*q(1,i)) 2*(q(3,i)*q(4,i)-q(2,i)*q(1,i)) q(1,i)^2-q(2,i)^2-q(3,i)^2+q(4,i)^2];
        euler(1,i-1) = atan2(C_dach{i-1}(2,3), C_dach{i-1}(3,3))/pi*180;
        euler(2,i-1) = asin(-C_dach{i-1}(1,3))/pi*180;
        euler(3,i-1) = atan2(C_dach{i-1}(1,2), C_dach{i-1}(1,1))/pi*180;
        [lon,lat,h(i)] = cart2ellip(pos(1,i),pos(2,i),pos(3,i),a,e);
        g_e = C_n2e(lat,lon) * g_n;
        v(:,i) = v(:,i-2) + (C_dach{i-2} * (3 * dv(:,i-1) - dv(:,i)) + 4 * C_dach{i-1} * (dv(:,i-1) + dv(:,i)) + C_dach{i} * (3 * dv(:,i) - dv(:,i-1))) / 6 ...
                 - (2 * Omega_ie * v(:,i-2) + Omega_ie * Omega_ie * pos(:,i-2) - g_e) * (2 * dt);
        pos(:,i+1) = pos(:,i) + v(:,i) * dt;
    end
end
figure
hold on
plot3(pos(1,:),pos(2,:),pos(3,:))
xlabel('x')
ylabel('y')
zlabel('z')
title('Trajektorie')
saveas(gcf,'trajektorie.png')

figure
hold on
subplot(3,1,1)
plot(t,pos(1,2:end))
xlabel('Zeit (s)')
ylabel('x (m)')
title('x')
subplot(3,1,2)
plot(t,pos(2,2:end))
xlabel('Zeit (s)')
ylabel('y (m)')
title('y')
subplot(3,1,3)
plot(t,pos(3,2:end))
xlabel('Zeit (s)')
ylabel('z (m)')
title('z')
saveas(gcf,'A2_pos.png')

figure
hold on
subplot(3,1,1)
plot(t,v(1,:))
xlabel('Zeit (s)')
ylabel('vx (/s)')
title('vx')
subplot(3,1,2)
plot(t,v(2,:))
xlabel('Zeit (s)')
ylabel('vy (/s)')
title('vy')
subplot(3,1,3)
plot(t,v(3,:))
xlabel('Zeit (s)')
ylabel('vz (/s)')
title('vz')
saveas(gcf,'A2_vel.png')

figure
hold on
subplot(3,1,1)
toplot = euler(1,:);
toplot(toplot<0) = toplot(toplot<0)+360;
plot(t(2:end),toplot)
xlabel('Zeit (s)')
ylabel('degree')
title('alpha')
subplot(3,1,2)
plot(t(2:end),euler(2,:))
xlabel('Zeit (s)')
ylabel('degree')
title('beta')
subplot(3,1,3)
plot(t(2:end),euler(3,:))
xlabel('Zeit (s)')
ylabel('degree')
title('gamma')
saveas(gcf,'A2_euler.png')
% 
% %% KF
% 
dt = 0.02;
% dt = 5;
% % % dt = 60;
dt = 120;
step = dt/(1/50);
len1 = floor(len/step);

tau = 10;
beta = 1 / tau;

pos_kf = zeros(3,len1);
v_kf = zeros(3,len1);


var_gauss_acc = 1e-5;
var_gauss_gyro = 1e-8;
x = zeros(15,1);
P = cell(len,1);
P{1} = diag([100;100;100;10;10;10;0.1;0.1;0.1;1;1;1;0.1;0.1;0.1]);
H = [eye(6),zeros(6,9)];
R = eye(6) * 1/1000;
code = 2;
if code == 1
    pos_kf = zeros(3,len1-1);
    v_kf = zeros(3,len1-1);
end
pos_kf(:,1) = pos(:,1);
v_kf(:,1) = v(:,1);
for i = 2:len1
    % F matrix
    F = zeros(15);
    F(1:3,4:6) = eye(3);
    F(4:6,1:3) = -Omega_ie * Omega_ie;
    F(4:6,4:6) = -2 * Omega_ie;
    F(4:6,7:9) = [0,-acc(3,(i-1)*step +1),acc(2,(i-1)*step +1); acc(3,(i-1)*step +1),...
        0,-acc(1,(i-1)*step +1); -acc(2,(i-1)*step +1),acc(1,(i-1)*step +1),0];
    [lon,lat,~] = cart2ellip(pos(1,(i-1)*step +1),pos(2,(i-1)*step +1),pos(3,(i-1)*step +1),a,e);
    C_p2e = C_n2e(lat,lon) * C_b2n(GyroX((i-1)*step +1),GyroY((i-1)*step +1),GyroZ((i-1)*step +1));
    F(4:6,10:12) = C_p2e;
    F(7:9,7:9) = -Omega_ie;
    F(7:9,13:15) = -C_p2e;
    F(10:12,10:12) = -eye(3) * beta;
    F(13:15,13:15) = -eye(3) * beta;
    % G matrix
    G = zeros(15,6);
    G(10:12,1:3) = eye(3) * sqrt(beta * var_gauss_acc);
    G(13:15,4:6) = eye(3) * sqrt(beta * var_gauss_gyro);
    A = [-F,G * G';zeros(15),F'] * dt;
    B = expm(A);
    Phi = B(16:30,16:30)';
    Q = Phi * B(1:15,16:30);
    
    xnnp = Phi * x;
    Pnnp = Phi * P{i-1} * Phi' + Q;
    K = Pnnp * H' * inv(H * Pnnp * H' + R);

    z1 = pos(:,1) - pos(:,i-1);
    z2 = [0;0;0] - v(:,i-1);
 
    z = [z1;z2];
    x_dach = xnnp + K * (z - H * xnnp);
    P{i} = (eye(15) - K * H) * Pnnp;

    v_kf(:,i) = v(:,(i-1)*step +1) + x_dach(4:6);  
    pos_kf(:,i) = pos(:,(i-1)*step +1) + x_dach(1:3);

end


figure
plot3(pos_kf(1,1:end-1)-pos(1,1),pos_kf(2,1:end-1)-pos(2,1),pos_kf(3,1:end-1)-pos(3,1));
%     hold on
%     plot3(pos(1,:),pos(2,:),pos(3,:))


figure
subplot(3,1,1)
plot(pos_kf(1,:)-pos(1,1:step:(i-1)*step +1))
title('pos diff -x')
subplot(3,1,2)
plot(pos_kf(2,:)-pos(2,1:step:(i-1)*step +1))
title('pos diff -y')
subplot(3,1,3)
plot(pos_kf(3,:)-pos(3,1:step:(i-1)*step +1))
title('pos diff -z')

figure
subplot(3,1,1)
plot(v_kf(1,:)-v(1,1:step:(i-1)*step +1))
title('vel diff -x')
subplot(3,1,2)
plot(v_kf(2,:)-v(2,1:step:(i-1)*step +1))
title('vel diff -y')
subplot(3,1,3)
plot(v_kf(3,:)-v(3,1:step:(i-1)*step +1))
title('vel diff -z')

figure
rr1 = sqrt(pos(1,:).^2 + pos(2,:).^2 + pos(3,:).^2);
rr2 = sqrt(pos_kf(1,:).^2 + pos_kf(2,:).^2 + pos_kf(3,:).^2);
subplot(2,1,1)
plot(rr1(1:end-1))
subplot(2,1,2)
plot(rr2(1:end-1))