% Die Funktion RK integriert numerisch Funktionen
% nach dem Runge-Kutta-Verfahren 3. Ordnung.

% Autor: Carsten Helfert 3318553
% Datum: 28.11.2020


function [q_2] = RK3_2(odefun, q_0, delta_t, omega_0, omega_1, omega_2)
    % RK3 Faktoren
    c1 = 1/3*delta_t;
    c2 = 4/3*delta_t;
    c3 = c1;
	
	% Berecknen der Elemente k
    k1 = odefun(q_0, omega_0);
    k2 = odefun(q_0 + k1*delta_t, omega_1);
    k3 = odefun(q_0 - k1*2*delta_t + k2 * 4*delta_t, omega_2);
	
	% Berechnen des Quaternions an der Stelle t_k+2
	q_2 = q_0 + c1*k1 + c2*k2 + c3*k3;
        
    % Normierung für Quaternionen
    q_2 = q_2 / sqrt(q_2'*q_2);
end

