function[X, Y, Z]=ellip2cart(lon,lat,h,a,e)
% lambda, phi, h can be vector
% if it is not given, a = 6378137 m
%     e = 0;
%     h = 0;

if nargin <5
    e = 0;
end 

if nargin < 4
    a = 6378137;
end

if nargin < 3
    h = 0;
end


% here we go
N = a ./ sqrt(1 - e^2 * sin(lat).^2);
X = (N + h) .* cos(lat) .* cos(lon);
Y = (N + h) .* cos(lat) .* sin(lon);
Z = (N .* (1 - e^2) + h) .* sin(lat);
end
