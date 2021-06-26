% Die Funktion odefun wird zur numerischen Integration
% nach dem Runge-Kutta-Verfahren verwendet.

% Autor: Carsten Helfert 3318553
% Datum: 28.11.2020


function [q_p] = odefun(q, omega)
    A = [0 omega(1) omega(2) omega(3)
         -omega(1) 0 omega(3) -omega(2)
         -omega(2) -omega(3) 0 omega(1)
         -omega(3) omega(2) -omega(1) 0];
	q_p = 1/2*A*q;
end

