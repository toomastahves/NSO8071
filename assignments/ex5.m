clear

% Load data
filename = 'global-reanalysis-phy-001-030-daily.nc';
%ncdisp(filename);
lon = ncread(filename,'longitude');
lat = ncread(filename,'latitude');
v = ncread(filename,'vo');
depth = ncread(filename,'depth');
time = ncread(filename,'time');
time = datetime(1950,1,1) + hours(time);

% a) Velocity across latitude 26N, longitude between -70 and -10
lat_index = find(lat == 26);
lon_index = find(lon >= -70 & lon <= -10);
z = v(lon_index, lat_index, :);
z(:,2) = []; % flatten matrix from 3D to 2D
z = rot90(z); % bottom down
x = double(lon(lon_index));
y = double(depth(1:49,:));

% Plot a)
h1 = pcolor(z);
set(h1, 'EdgeColor', 'none');
title('Northward velocity vs depth');

xlabel('Longitude');
x_label = x(1:72:end);
c = ismember(x, x_label);
x_index = find(c);
set(gca,'xtick', x_index);
set(gca,'xticklabel', x_label);

ylabel('Depth (m)');
y_label = y(1:5:end);
c = ismember(y, y_label);
y_index = find(c);
set(gca,'ytick', y_index);
set(gca,'yticklabel', flip(round(y_label + 2000)));

yyaxis right;
ylabel('Northward sea water velocity (m/s)');
set(gca,'ytick', []);

saveas(h1, 'ex5_velocity','jpg');

% b) Integrate v-velocity over latitudes 30-40N
lat_index = find(lat >= 30 & lat <= 40);
z = v(lon_index, lat_index, :) * (1852*60*mean(diff(lat(lat_index))));
z_sum = sum(z,2); % integrate over latitudes
z_sum(:,2) = []; % flatten matrix from 3D to 2D
z_sum = rot90(z_sum); % bottom down
x = double(lon(lon_index));
y = double(depth(1:49,:));

% Plot b)
h2 = pcolor(z_sum);
set(h2, 'EdgeColor', 'none');
title('Northward acceleration vs depth');

xlabel('Longitude');
x_label = x(1:72:end);
c = ismember(x, x_label);
x_index = find(c);
set(gca,'xtick', x_index);
set(gca,'xticklabel', x_label);

ylabel('Depth (m)');
y_label = y(1:5:end);
c = ismember(y, y_label);
y_index = find(c);
set(gca,'ytick', y_index);
set(gca,'yticklabel', flip(round(y_label + 2000)));

yyaxis right;
ylabel('Northward sea water acceleration (m/s^2)');
set(gca,'ytick', []);

saveas(h2, 'ex5_acceleration','jpg');
