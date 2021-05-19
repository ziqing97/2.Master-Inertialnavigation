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
w_epp = w_ipp;     
q0 = [1 0 0 0];
q = zeros(length(t),4);
q = [q0;q];
dt = 1;
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
DCM = DCM(2:end);

a_e = zeros(3,2000);
v_e = zeros(length(t),3)';
x_e = zeros(length(t),3)';
for i = 1 : 2000
    if norm(x_e(:,i)) == 0
        g = [0;0;0];
    else
        g = -9.81 * x_e(:,i) / norm(x_e(:,i));
    end
    a_e(:,i) = DCM{i} * a(:,i) + g;
    if i > 1
        v_e(:,i) = v_e(:,i-1) + a_e(:,i-1) * dt;
    end
    if i > 1
        x_e(:,i+1) = x_e(:,i) + v_e(:,i-1) *dt;
    end
end

r_e = sqrt(x_e(1,:).^2 + x_e(2,:).^2 + x_e(3,:).^2);
index = find(r_e > 6378000);

figure;
plot3(x_e(1,:),x_e(2,:),x_e(3,:),'LineWidth',4);
hold on
scatter3(0,0,0)
title('Trajektorie')
grid on
% saveas(gcf,'Trajektorie.png')

v_start = v_e(:,end);
x_start = x_e(:,end);
x_second = x_e(:,end-1);


%% Aufgabe 2
a = flip(a,2);
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

DCM = DCM(2:end);

v_e = cell(1,8);
x_e = cell(1,8);
a_e = cell(1,8);


for j=1:8
    v_e{j} = zeros(length(t),3)';
    x_e{j} = zeros(length(t),3)';
    a_e{j} = zeros(3,2000);
    an(1,:) = a(1,:);
    an(2,:) = a(2,:); 
    an(3,:) = a(3,:); 
    an(1,:) = a(1,:) + bias_a(1,j);
    an(2,:) = a(2,:) + bias_a(2,j);
    an(3,:) = a(3,:) + bias_a(3,j);
    v_e{j}(:,1) = v_start;
    x_e{j}(:,1) = x_start;
    x_e{j}(:,2) = x_second;
    for i = 1 : 2000
        g = -9.81 * x_e{j}(:,i) / norm(x_e{j}(:,i));
        a_e{j}(:,i) = DCM{i} * an(:,i) + g;
        if i > 1
            v_e{j}(:,i) = v_e{j}(:,i-1) - a_e{j}(:,i-1) * dt;
        end
        if i > 1
            x_e{j}(:,i+1) = x_e{j}(:,i) - v_e{j}(:,i-1) *dt; 
        end
    end
    figure
    plot3(x_e{j}(1,:),x_e{j}(2,:),x_e{j}(3,:),'LineWidth',4);
    hold on
    grid on
    title(bias_a(:,j)')
%     saveas(gcf,num2str(j),'png')
end

