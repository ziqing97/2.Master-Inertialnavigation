function[C] = C_n2e(lat,lon)
    % lat lon in grad
    C = [-sin(lat) * cos(lon), -sin(lon), -cos(lat) * cos(lon);
         -sin(lat) * sin(lon), cos(lon), -cos(lat) * sin(lon);
         cos(lat), 0, -sin(lat)];
end