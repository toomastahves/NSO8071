clear

x = (-10 : 0.05: 10);
y = (1./(2*x)).*cos(2*x) + exp(-2*x).* sin(4*x);

fig = figure(1);
plot(x,y);
title('plot');
xlabel('x');
ylabel('f(x)');
legend({'f(x)'});
grid on;
saveas(fig, 'ex0_Figure1_ToomasTahves.png');

x = linspace(0,500);
y = linspace(-100,100);
k = 2 * pi / 100;
fi0 = 25 * k;
l = 2 * pi / 50;

[X,Y] = meshgrid(x,y);
Z = 2 * cos(k * X + fi0) .* cos(l * Y);
v = [-1,0,1];

[M,c] = contourf(X,Y,Z,v);
c.LineColor = 'black';
colorbar;
xlabel('X');
ylabel('Y');

f = gcf;
exportgraphics(f,'ex0_Figure2_ToomasTahves.png','Resolution',300);
