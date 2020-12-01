clear

% Define initial variables
g = 9.81; % m/s^2
alpha = 1e-5; % rad
f = 1.45*1e-4; % 1/s
x0 = 0; % m
y0 = 0; % m
u0 = 0; % m/s
v0 = 0; % m/s
dt = 360; % s
t = (0:dt:24*3600); % s
len = length(t);

% Analytical solution
A = g * sin(alpha) / f;
u_an = A * sin(f*t);
v_an = -A + A*cos(f*t);
x_an = -A / f*cos(f*t);
y_an = -A * t + A / f*sin(f*t);
velocity_an = sqrt(u_an.^2 + v_an.^2);
trajectory_an = sqrt(x_an.^2 + y_an.^2);

% Explicit solution
x_expl = zeros(1,len);
y_expl = zeros(1,len);
u_expl = zeros(1,len);
v_expl = zeros(1,len);
for n = 1:len-1
    u_expl(n+1) = u_expl(n) + dt * (f*v_expl(n) + g*sin(alpha));
    v_expl(n+1) = v_expl(n) + dt * (-f*u_expl(n));
    x_expl(n+1) = x_expl(n) + dt * u_expl(n);
    y_expl(n+1) = y_expl(n) + dt * v_expl(n);
end
velocity_expl = sqrt(u_expl.^2 + v_expl.^2);
trajectory_expl = sqrt(x_expl.^2 + y_expl.^2);

% Implicit solution
x_impl = zeros(1,len);
y_impl = zeros(1,len);
u_impl = zeros(1,len);
v_impl = zeros(1,len);
for n = 1:len-1
    coeff = 1 + dt^2 * f^2;
    u_impl(n+1) = (u_impl(n) + dt * (f*v_impl(n) + g*sin(alpha))) / coeff;
    v_impl(n+1) = (v_impl(n) + dt * (-f*u_impl(n+1)));
    x_impl(n+1) = x_impl(n) + dt * u_impl(n+1);
    y_impl(n+1) = y_impl(n) + dt * v_impl(n+1);
end
velocity_impl = sqrt(u_impl.^2 + v_impl.^2);
trajectory_impl = sqrt(x_impl.^2 + y_impl.^2);

% Plotting
h1 = figure;
hold on;
grid on;
title('Trajectory vs time');
xlabel('Time') 
ylabel('Trajectory');
plot(t, trajectory_expl, 'DisplayName','Explicit');
plot(t, trajectory_impl, 'DisplayName','Implicit');
plot(t, trajectory_an, 'DisplayName','Analytical');
legend('Explicit','Implicit','Analytical');
hold off;
saveas(h1, 'ex1_trajectory','jpg');

h2 = figure;
hold on;
grid on;
title('Velocity vs time');
xlabel('Time') 
ylabel('Velocity');
plot(t, velocity_expl, 'DisplayName','Explicit');
plot(t, velocity_impl, 'DisplayName','Implicit');
plot(t, velocity_an, 'DisplayName','Analytical');
legend('Explicit','Implicit','Analytical');
hold off;
saveas(h2, 'ex1_velocity','jpg');

% Error calculation
err_impl_vel = velocity_impl - velocity_an;
err_expl_vel = velocity_expl - velocity_an;
err_impl_traj = trajectory_impl - velocity_an;
err_expl_traj = trajectory_expl - velocity_an;

% Error plotting
h3 = figure;
hold on;
grid on;
title('Velocity error for numerical vs analytical solution');
xlabel('Time') 
ylabel('Velocity error vs time');
plot(t, err_expl_vel, 'DisplayName','Explicit');
plot(t, err_impl_vel, 'DisplayName','Implicit');
legend('Explicit','Implicit');
hold off;
saveas(h3, 'ex1_velocity_error','jpg');

h4 = figure;
hold on;
grid on;
title('Trajectory error for numerical vs analytical solution');
xlabel('Time') 
ylabel('Trajectory error');
plot(t, err_expl_traj, 'DisplayName','Explicit');
plot(t, err_impl_traj, 'DisplayName','Implicit');
legend('Explicit','Implicit');
hold off;
saveas(h4, 'ex1_trajectory_error','jpg');
