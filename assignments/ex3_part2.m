clear;

% Load data
load('ex3_data.mat', 'bath', 'distance');

% Grid size
M = 72;
N = 3;

% PARAMETERS
h1 = figure;
f = 10^(-4); % Coriolis parameter
g = 10;  % gravity
h0 = 67.3; % reference depth
R = sqrt(g*h0)/f; % deformation radius

Lx = 1200*R; % Domain size
Ly = 1000*R; % Domain size
%dx = Lx/(M-2); % Periodic conditions on x require a ghost point on each side
%dy = Ly/(N-2); % Impermeability in y direction two points are used for land
dx = 2300; % Periodic conditions on x require a ghost point on each side
dy = 2300; % Impermeability in y direction two points are used for land

%TF = 1*Lx/sqrt(g*h0); % Final time here one displacement 
TF = 5000;
dt = 0.1/sqrt(g*h0)/sqrt((1/dx^2+1/dy^2)); % Time step based on gravity wave speed
t = 0;
Nsteps = floor(TF/dt); % Number of time steps as a function of time step and Final time

% Coordinates of cell centers
for i=1:M
    x(i) = (i-1/2)*dx;
end
% Zero in j=1+1/2. First line is a land point
for j=1:N
    y(j) = (j-3/2)*dy;
end

% RESOLUTION AND PLOT
% Velocities
u = zeros(M,N);
v = zeros(M,N);
% New velocities
unp=u;
vnp=v;
% Transports
hu=u;
hv=v;
% Elevation
zeta = zeros(M,N);
% Topography
h = h0*ones(M,N);
% Land-sea mask
issea = ones(M,N);

% Staggered C-grid:
% u(i,j) is in i-1/2
% v(i,j) is in j-1/2

kx = 2*pi/Lx; % Initial situation: a full wave
% Small amplitude. For larger amplitudes, nonlinearities appear
% Original code replaced with data calculated in ex3_part1.m
for j=1:N
    for i=1:M
        h(i,j) = -bath(i);            
    end
end

for i=1:M
    for j=1:N
        zeta(i,j) = 0; 
    end
end

htot=h+zeta; % Total water depth including elevation
% Solid boundary in south and north
issea(:,1) = 0;
issea(:,N) = 0;
issea(h<0) = 0;
zeta(1:5,:) = 0.1;

% Show initial elevation
%figure(1);
%pcolor(x,y,h');
%shading flat;
%title ('To start press any key');
%colormap;
%pause;

% Make Nsteps time steps
seaLevels = zeros(Nsteps,1); % Sealevels near Ecuador
for nn=1:Nsteps
    % Defining transports from velocities using actual total depth (nonlinear)
    for i=2:M
        for j=1:N
            hu(i,j) = u(i,j)*(htot(i,j)+htot(i-1,j))/2;
        end
    end
    hu(1,:)=0;
    hu(M,:)=0;
    for i=1:M
        for j=2:N
            hv(i,j) = v(i,j)*(htot(i,j)+htot(i,j-1))/2;
        end
    end
    hv(:,1)=0;

    % New elevation, updated immediately since local and never needing any neighbour zeta
    for i=2:M-1
        for j=2:N-1
            % divergence of transport leads to elevation change
            zeta(i,j) = zeta(i,j) - dt/dx * (hu(i+1,j)-hu(i,j)) - dt/dy * (hv(i,j+1)-hv(i,j));
        end
    end

    % Periodicity
    %zeta(1,:)=zeta(M-1,:);
    %zeta(M,:)=zeta(2,:);
    zeta(1,:)=0;
    zeta(M,:)=0;
    zeta(:,1)=0;
    zeta(:,N)=0;
    % Total height
    htot = h+zeta;
    % New velocity values with fractional-step approach on coriolis
    if mod(nn,2) == 1
        multipas = [1 2];
    elseif mod(nn,2) == 0
        multipas = [2 1];
    end

    % Loop
    for k=1:2
        if multipas(k) == 1
            for i=2:M-1
                for j=1:N-1
                    % Coriolis
                    work = u(i,j) + f*dt/4 * (v(i,j+1) + v(i,j) + v(i-1,j+1) + v(i-1,j));
                    % Elevation gradientz
                    unp(i,j) = work - g*dt/dx * (zeta(i,j) - zeta(i-1,j));
                    % Applying a mask on sea-land interfaces
                    unp(i,j)=unp(i,j)*issea(i,j)*issea(i-1,j);
                end
            end
            % Update and periodic conditions
            u = unp;
            %u(1,:)=u(M-1,:);
            %u(M,:)=u(2,:);
            u(1,:)=0; 
            u(M,:)=0;
            u(:,N)=0;
        elseif multipas(k) == 2
            for i=2:M-1
                for j=2:N
                    % Coriolis
                    work = v(i,j) - f*dt/4 * (u(i+1,j) + u(i+1,j-1) + u(i,j) + u(i,j-1));
                    % and elevation gradient
                    vnp(i,j) = work - g*dt/dy * (zeta(i,j) - zeta(i,j-1));
                    % applying a mask on sea-land interfaces
                    vnp(i,j)=vnp(i,j)*issea(i,j)*issea(i,j-1);
                end
            end
            % update and periodic conditions
            v = vnp;
            %v(1,:)=v(M-1,:);
            %v(M,:)=v(2,:);
            v(1,:)=0;
            v(:,1)=0;
            v(M,:)=0;
        end
    end

    seaLevels(nn) = zeta(71,2);
    t=t+dt; % update time
    % Every once in a while make a plot
    if mod(nn,20)==0
        plot(x, zeta(:,2));
        title(['Elevation at time: ',num2str(t)]);
        xlabel('Distance (m)');
        ylabel('Wave height (m)');
        axis([1000 165000 -0.1 0.15])
        grid on;
        drawnow;
        frame = getframe(h1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if nn == 20
            imwrite(imind,cm,'ex3_animation.gif','gif', 'Loopcount',inf); 
        else
            imwrite(imind,cm,'ex3_animation.gif','gif','WriteMode','append'); 
        end
        pause(0.01);
    end
end

% Plot sea levels near Ecuador
steps = 1:Nsteps;
h2 = figure;
plot(steps*dt/3600, seaLevels);
xlabel('Time (h)');
ylabel('Sea surface elevation (m)');
title('Sea elevation near Ecuador');
grid on;
saveas(h2, 'ex3_sealevels','png');
