function[C] = C_b2n(yaw,pitch,roll)
    Y = yaw;
    P = pitch;
    R = roll;
    C = [cos(Y) * cos(P), cos(Y) * sin(P) * sin(R) - sin(Y) * cos(R), cos(Y) * sin(P) * cos(R) + sin(Y) * sin(R);
        sin(Y) * cos(P), sin(Y) * sin(P) * sin(R) + cos(Y) * cos(R), sin(Y) * sin(P) * cos(R) - cos(Y) * sin(R);
        -sin(P), cos(P) * sin(R), cos(P) * cos(R)];
end