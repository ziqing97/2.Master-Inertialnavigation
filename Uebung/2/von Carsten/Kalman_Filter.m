% X_d - Vorherige ausgeglichene Beobachtung
% Phi - Übergangsmatrix
% P - Kovarianz des Zustandes
% Q - Matrix mit Prozessrauschen
% H - Verbindung zwischen Zustand und Beobachtung
% R - Matrix mit Messrauschen

function [x, P] = Kalman_Filter(x_0, pos_error, var_pos, vel_error, ...
                                var_vel, P, var_acc, var_gyro, ...
                                acc, DCM, Omega_ie, tau, delta_t)
%% Initialisierungen
    % F-Matrix
    beta = 1/tau;
    
    A_e = [0 -acc(3) acc(2)
           acc(3) 0 -acc(1)
           -acc(2) acc(1) 0];
    Fa = -eye(3)*beta;
    Fw = Fa;

    F = [zeros(3,3) eye(3) zeros(3,3) zeros(3,3) zeros(3,3)
         -Omega_ie*Omega_ie -2*Omega_ie A_e DCM zeros(3,3)
         zeros(3,3) zeros(3,3) -Omega_ie zeros(3,3) -DCM
         zeros(3,3) zeros(3,3) zeros(3,3) Fa zeros(3,3)
         zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) Fw];
    
    Ga = [sqrt(beta*var_acc) 0 0
          0 sqrt(beta*var_acc) 0
          0 0 sqrt(beta*var_acc)];
    Gw = [sqrt(beta*var_gyro) 0 0
          0 sqrt(beta*var_gyro) 0
          0 0 sqrt(beta*var_gyro)];
    
    G = [zeros(3,3), zeros(3,3); 
         zeros(3,3), zeros(3,3); 
         zeros(3,3), zeros(3,3);
         Ga, zeros(3,3);
         zeros(3,3), Gw];
    
   GWG = G*G';
   
   % Kochbuchrezept
   A = [-F GWG
        zeros(15,15) F']*delta_t;
   B = expm(A);
   Phi = B(16:end, 16:end)';
   Q = Phi * B(1:15, 16:end);
   
   H = [eye(6) zeros(6,9)];
   R = eye(6).*[var_pos var_pos var_pos var_vel var_vel var_vel];
   z_n = [pos_error
          vel_error];

%% Kalman Filter
    % X_d_n|n-1
    x = Phi*x_0;
    
    % P_n|n-1
    P = Phi*P*Phi'+Q;
    
    % K_n
    K = P*H'*inv(H*P*H'+R);
    
    % Update mit z_n
    % X_d_n|n
    x = x+K*(z_n-H*x);
    
    % Kovarianz P_n|n
    P = (eye(length(H))-K*H)*P;
end

