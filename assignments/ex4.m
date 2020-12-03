clear

% For plotting
h = figure;
hold on;
grid on;
title('Earth temperature vs albedo');
xlabel('Temperature (K)');
ylabel('Albedo');

% Variables
alpha = 0; % Albedo - describes how much sunlight Earth reflects.
sigma = 5.67 * 1e-8; % Stefan-Boltzmann constant - total intensity.
I = 344; % Incident - short-wave radiation from the sun.
% Base function - balance between amount of emitted/absorbed energy:
% (1 - alpha) * I  = sigma * T^4

% Case 1:
% When temperature on Earth is T1 = 235 K,
% then alpha is 0.5 (high albedo - there is ice and clouds reflecting sun).
alpha1 = 0.5;
T1_mean = nthroot((1-alpha1)*I/sigma, 4); % 234 K
line([200 250],[alpha1 alpha1],'LineWidth',1);

% Case 2:
% When temperature is between 250 and 270 K, then albedo has functional
% dependency of (270 - T) / 40. After solving balance equation and summing
% temperatures, then mean is 259 K.
T2 = 250:1:270;
alpha2 = (270 - T2) / 40;
syms T; % Define symbol for solving equation
eqnLeft = T; % Left side of equation
eqnRight = nthroot((1-(270 - T) / 40)*I/sigma, 4); % Right side of equation
initialGuess = 300; % Initial guess where solution converges to 260 K.
T_mean = vpasolve(eqnLeft == eqnRight, T, initialGuess); % 260 K
plot(T2, alpha2);

% Case 3:
% When temperature on Earth is T3 = 280 K,
% then alpha is 0 (low albedo - no ice and clouds reflecting sun).
alpha3 = 0;
T3_mean = nthroot((1-alpha3)*I/sigma, 4); % 279 K
line([270 300],[alpha3 alpha3],'LineWidth',1);

% Save plot to file
saveas(h, 'ex4','jpg');

% Summary: 
% In first case average temperature is 234 K, 
% in second case 259 K 
% and in third case 279 K.
