clear all
lat = [20:1:30];
lon = [80:1:90];

% pre-initialize
lat_mat = zeros(11,11);
lon_mat = zeros(11,11);

for ii = 1:length(lat)
    for jj = 1:length(lon)
        lat_mat(ii,jj) = lat(ii);
        lon_mat(ii,jj) = lon(jj);
    %fprintf(' lat = %.2f lon = %.2f \n', lat(ii),lon(jj))
    end 
end


