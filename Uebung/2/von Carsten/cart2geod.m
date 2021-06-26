function [lat, lon, alt] = cart2geod(pos)
a = 6378137;
f = 1.0/298.257223563;
e2 = (2-f)*f;
lon = atan2(pos(2),pos(1));
lat = atan2(pos(3)/sqrt(pos(1)^2+pos(2)^2),0.01);
N = a/sqrt(1-e2*sin(lat).^2);
alt = sqrt(pos(1)^2+pos(2)^2)/cos(lat)-N;
alt_vgl = 0;
while abs(alt-alt_vgl) > 1e-6
	  alt_vgl  = alt;
	  lat = atan2(pos(3)/sqrt(pos(1)^2+pos(2)^2),1-e2*N/(N+alt));
	  N = a/sqrt(1-e2*sin(lat)^2);
	  alt = sqrt(pos(1)^2+pos(2)^2)/cos(lat)-N;
end
end