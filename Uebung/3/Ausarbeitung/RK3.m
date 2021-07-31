function[qk] = RK3(qk_m2,omega_m2, omega_m1, omegak,dt)
    k1 = 0.5 * A(omega_m2) * qk_m2;
    k2 = 0.5 * A(omega_m1) * (qk_m2 + k1 * dt);
    k3 = 0.5 * A(omegak) * (qk_m2 - k1 * 2 * dt + k2 * 4 * dt);
    qk = qk_m2 + dt/3 * (k1 + 4 * k2 + k3);
end