clear;

% Load daily data
data = 'reanalysis_dailymeans_20180626.nc';
longitude = ncread(data,'longitude');
latitude = ncread(data,'latitude');
depth = ncread(data, 'depth');

% Get latitude indexes
start_index_lon = find(longitude > 23.13,1,'first');
start_index_lat = find(latitude > 59.07,1,'first');
latitude_range = latitude(start_index_lat:start_index_lat+22);

% Select data corresponding to ranges
start = [start_index_lon, start_index_lat, 1, 1];
count = [1, 23, 56, 1];
temperature = squeeze(ncread(data, 'thetao', start, count));
salt = squeeze(ncread(data, 'so' , start, count));

% Select depth range
depth_indexes = squeeze(find(depth <= 90));
depth = depth(depth_indexes);

% Calculate density
density = 1028*(1-0.00017*(temperature(:,depth_indexes)-10)+ 0.00076*(salt(:,depth_indexes)-35))';

% Plot
h1 = pcolor(latitude_range, depth, density);
title('Density distribution along transect for 26. June');
xlabel('Latitude');
ylabel('Depth (m)');
c = colorbar;
c.Label.String = 'Density (kg/m^3)';
set(h1, 'EdgeColor', 'none');
set(gca,'YDir','reverse')
saveas(h1, 'ex7','jpg');

% Just calculations without plot
% Derivative of density
drho = diff(density(depth_indexes,:),1);
% Estimate depth of thermocline h0
h0 = depth(squeeze(find(drho(:,floor(length(latitude_range)/2)) == nanmax(drho(:,floor(length(latitude_range)/2))))));
% Upper layer
rho_upper = nanmean(density(squeeze(find(depth < h0)),:),'all');
% Lower layer
rho_lower = nanmean(density(squeeze(find(depth > h0)),:),'all');
% Calculate reduced gravity
reduced_gravity = 9.81*(rho_lower - rho_upper) / 1028;
