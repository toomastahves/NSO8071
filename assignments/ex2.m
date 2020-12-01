clear

% ncdisp('data.nc');
%lon = ncread('data.nc','longitude'); % 2
%lat = ncread('data.nc','latitude'); % 4
u = ncread('data.nc','u10');
v = ncread('data.nc','v10');
pressure = ncread('data.nc','msl');
time = ncread('data.nc','time');
time = datetime(1900,1,1) + hours(time);
temp = ncread('data.nc','t2m') - 273.15;
temp = temp(2,4,:); % select values at 22.25 lon and 59.25 lat
temp = temp(1,:)';
h = figure;
hold on;
grid on;
title('Temperature vs time');
xlabel('Time (Months)') 
ylabel('Temperature (Celsius)');
plot(time, temp);
saveas(h, 'ex2_temperature','jpg');

lon = 2; % 22.25
lat = 4; % 59.25
f = 4*(pi/(24*60*60))*sin(59.25*pi/180);
dy = 0.25 * 60 * 1852;
dx = 0.25 * 60 * 1852 * cos(59.25*pi/180);
dt = 60 * 60; % h
timespan = (7993:1:8016); % 30.11

len = length(timespan);
acc_local = zeros(1,len);
acc_coriolis = zeros(1,len);
acc_pressure = zeros(1,len);
acc_advection = zeros(1,len);
acc_other = zeros(1,len);

for n = 2:len-1
    acc_coriolis(n) = -1*f*v(lon, lat, timespan(n));
    acc_local(n) = (u(lon, lat, timespan(n+1)) - u(lon, lat, timespan(n-1))) / (2*dt);
    acc_pressure(n) = -1*(pressure(lon+1, lat, timespan(n))-pressure(lon-1,lat,timespan(n)))/(2*dx);
    
    adv_u = u(lon,lat,timespan(n)) * (u(lon+1,lat ,timespan(n)) - u(lon-1,lat,timespan(n))) / (2*dx);
    adv_v = v(lon,lat,timespan(n)) * (v(lon, lat,timespan(n)) - v(lon+1,lat,timespan(n))) / (2*dy);
    acc_advection(n) = adv_u + adv_v;
    
    acc_other = acc_coriolis + acc_local + acc_advection - acc_pressure;
end

h1 = figure;
hold on;
grid on;
title('Acceleration vs time');
xlabel('Time (hours)') 
ylabel('Acceleration (m/s^2)');
time_plot = time(timespan);
time_plot = time_plot(2:end-1);
plot(time_plot, acc_coriolis(2:end-1), 'DisplayName','Coriolis');
plot(time_plot, acc_local(2:end-1), 'DisplayName','Local');
plot(time_plot, acc_pressure(2:end-1), 'DisplayName','Pressure');
plot(time_plot, acc_advection(2:end-1), 'DisplayName','Advection');
plot(time_plot, acc_other(2:end-1), 'DisplayName','Other');
legend('Coriolis','Local','Pressure','Advection','Other');
hold off;
saveas(h1, 'ex2_acceleration','jpg');
