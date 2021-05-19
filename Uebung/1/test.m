% INav Ue1
% Ziqing Yu
% 3218051
%% Initial
clc
close all
clearvars

%% data
imudata = importfile('imu-data.txt', 2, 2001);
t = imudata(:,1);
a = imudata(:,2:4);
a = a';
w_ipp = imudata(:,5:7); % delta alpha

%% Aufgabe 1
w_epp = w_ipp;     % w_e = 0, delta_beta = delta_alpha
q0 = [1 0 0 0];
q = zeros(length(t),4);
q = [q0;q0;q];
dt = 1;
q = q';
for i = 1:length(t)
    A = [0, w_epp(i,1), w_epp(i,2), w_epp(i,3);
         -w_epp(i,1), 0, w_epp(i,3), -w_epp(i,2);
         -w_epp(i,2), -w_epp(i,3), 0, w_epp(i,1);
         -w_epp(i,3), w_epp(i,2), -w_epp(i,1), 0];
    q(:,i+2) = expm(0.5 * A * dt) * q(:,i+1);
end

% DCM
DCM = cell(2002,1);
for i = 1:length(t)+2
    q0 = q(1,i);
    q1 = q(2,i);
    q2 = q(3,i);
    q3 = q(4,i);
    DCM{i} = [q0^2 + q1^2 - q2^2 - q3^2, 2 * (q1 * q2 + q3 * q0), 2 * (q1 * q3 - q2 * q0);
               2 * (q1 * q2 - q3 * q0), q0^2 - q1^2 + q2^2 - q3^2, 2 * (q2 * q3 + q1 * q0);
               2 * (q1 * q3 + q2 * q0), 2 * (q2 * q3 - q1 * q0), q0^2 - q1^2 - q2^2 + q3^2];
end


% Position
v_e = zeros(length(t)+1,3)';
x_e = zeros(length(t),3)';
r_e = zeros(length(t),1);

for i = 3:length(t)+1
    if norm(x_e(:,i-2)) == 0
        g = 0;
    else
        g = -9.81 * x_e(:,i-2) / norm(x_e(:,i-2)); % g
    end
    v_e(:,i) = v_e(:,i-2) + ...
             (DCM{i-1} * (3 * a(:,i-2) - a(:,i-1)) + ...
             4 * DCM{i} * (a(:,i-2) + a(:,i-1)) + ...
             DCM{i+1} * (3 * a(:,i-1) - a(:,i-2))) / 6 - ...
             (-g) * 2 * dt;
         
     x_e(:,i-1) = x_e(:,i-2) + v_e(:,i) * dt;
     r_e(i-1) = norm(x_e(:,i-1));
end

figure;
plot3(x_e(1,:),x_e(2,:),x_e(3,:),'LineWidth',4);
title('Trajektorie')
grid on


%% Aufgabe 2
a = -flip(a,2);
w_epp = -flip(w_epp,1);

bias_a = 0.1 * [1,1,1;
                1,1,-1;
                1,-1,1;
                1,-1,-1;
                -1,1,1;
                -1,1,-1;
                -1,-1,1;
                -1,-1,-1]'; % m/s^2


q = zeros(length(t),4);
dt = 1;
q0 = [1 0 0 0];
q = [q0;q];
q = q';
for i = 1:length(t)
    A = [0, w_epp(i,1), w_epp(i,2), w_epp(i,3);
         -w_epp(i,1), 0, w_epp(i,3), -w_epp(i,2);
         -w_epp(i,2), -w_epp(i,3), 0, w_epp(i,1);
         -w_epp(i,3), w_epp(i,2), -w_epp(i,1), 0];
    q(:,i+1) = expm(0.5 * A * dt) * q(:,i);
end

% DCM
DCM = cell(2001,1);
for i = 1:length(t)+1
    q0 = q(1,i);
    q1 = q(2,i);
    q2 = q(3,i);
    q3 = q(4,i);
    DCM{i} = [q0^2 + q1^2 - q2^2 - q3^2, 2 * (q1 * q2 + q3 * q0), 2 * (q1 * q3 - q2 * q0);
               2 * (q1 * q2 - q3 * q0), q0^2 - q1^2 + q2^2 - q3^2, 2 * (q2 * q3 + q1 * q0);
               2 * (q1 * q3 + q2 * q0), 2 * (q2 * q3 - q1 * q0), q0^2 - q1^2 - q2^2 + q3^2];
end

v_start = v_e(:,end);
x_start = x_e(:,end-1:end);

v_e = cell(1,8);
x_e = cell(1,8);
str = zeros(3,8);
for i = 1:8
    a(1,:) = a(1,:) + bias_a(1,i);
    a(2,:) = a(2,:) + bias_a(2,i);
    a(3,:) = a(3,:) + bias_a(3,i);
    % Position
    v_e{i} = zeros(length(t),3)';
    x_e{i} = zeros(length(t),3)';
    
    v_e{i}(:,1) = [0;0;0];
    v_e{i}(:,2) = -v_start(:,1);
    x_e{i}(:,1) = x_start(:,2);
    x_e{i}(:,2) = x_start(:,1);

    x_e{i}(:,3) = x_e{i}(:,2) + v_e{i}(:,2) * dt;
    for j = 4:length(t)
        g = -9.81 * x_e{i}(:,j-1) / norm(x_e{i}(:,j-1)); % g
        v_e{i}(:,j) = v_e{i}(:,j-2) + ...
                 (DCM{j-2} * (3 * a(:,j-1) - a(:,j)) + ...
                 4 * DCM{j-1} * (a(:,j-1) + a(:,j)) + ...
                 DCM{j} * (3 * a(:,j) - a(:,j-1))) / 6 - ...
                 (-g) * 2 * dt;

         x_e{i}(:,j) = x_e{i}(:,j-1) + v_e{i}(:,j-1) * dt;
    end
    str(:,i) = x_e{i}(:,end);
    figure
    plot3(x_e{i}(1,:),x_e{i}(2,:),x_e{i}(3,:),'LineWidth',4);
    title(bias_a(:,i)')
    grid on

end




