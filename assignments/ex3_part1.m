clear;

% Load data from file
lat = ncread('ETOPO2v2c_f4.nc','y');
lon = ncread('ETOPO2v2c_f4.nc','x');
z = ncread('ETOPO2v2c_f4.nc','z');

% Extract sub-region
figure(1);
lat_indexes = 2400:3100;
lon_indexes = 2400:3100;
zsub = z(lon_indexes, lat_indexes);
lats = lat(lat_indexes);
lons = lon(lon_indexes);

% Plot map
pcolor(zsub);
hold on;
contour(zsub,[0 0],'edgecolor','k');
shading flat;
xlabel('Longitude');
ylabel('Latitude');
title('Isla Isabela and coast of Ecuador');
text(278, 274, 'Isla Isabela');
text(207, 595, 'Ecuador');

% Plot line
xind = 278:-1:207;
yind=round(linspace(274, 595, length(xind)));
h1 = plot(xind, yind, '-k');
saveas(h1, 'ex3_map','jpg');

% Calcualte bathymetry
bath = xind*0;
distance = bath*0;
for i=1:length(xind)
   bath(i)=zsub(yind(i), xind(i));
   if i>1
       dx=(lons(xind(i))-lons(xind(i-1))) *cosd(lats(yind(i)))*60*1851;
       dy=(lats(yind(i))-lats(yind(i-1)))*60*1851;
       distance(i)=distance(i-1)+ abs( dx+1i*dy );
   end
end

% Plot bathymetry
figure(2);
h2 = plot(distance*1e-3,bath);
title('Bathymetry between Isla Isabela and coast of Ecuador');
xlabel('Distance (km)');
ylabel('Depth (m)');
grid on;
saveas(h2, 'ex3_bathymetry','jpg');

% Save data
save('ex3_data.mat', 'bath', 'distance');
