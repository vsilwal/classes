clc
clear all

% variable input method
input_method = 1;

if input_method == 1
    % Initialization
    lat = [20:1:33];
    lon = [80:1:88];
else
    % ALTERNATE: USER INPUT
    min_lat = input('Minimum latitute: ');
    max_lat = input('Maximum latitute: ');
    min_lon = input('Minimum longitude: ');
    max_lon = input('Maximum longitude: ');
    dleta_lat = 1;
    delta_lon = 1;
    lat = min_lat:delta_lat:max_lat;
    lon = min_lon:delta_lon:max_lon;
end

% pre-initialize
lat_mat = nan(length(lat),length(lon));
lon_mat = nan(length(lat),length(lon));

[m,n] = size(lon_mat);

for ii = 1:m
    for jj = 1:n
        lat_mat(ii,jj) = lat(ii);
        lon_mat(ii,jj) = lon(jj);
        fprintf('lat = %d lon = %d \n', lat(ii),lon(jj))
    end 
end

