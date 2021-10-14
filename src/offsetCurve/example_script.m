%% EXAMPLE SCRIPT

% Define an x-y curve
x = [0 3 2 4 4 5 6 5 7];
y = [8 12 20 22 14 8 16 20 20];

plot(x, y, '-k');
hold on;

% Compute offset curve by 0.3
[xo, yo] = offsetCurve(x, y, 0.3);
plot(xo, yo, '-b');

% Compute offset curve by -0.8
[xo2, yo2] = offsetCurve(x, y, -0.8);
plot(xo2, yo2, '-r');

% Compute offset curve by 6 points respecting proportions of current axis
drawnow();
set(gca, 'XLimMode', 'manual'); % First lock the X and Y lims
set(gca, 'YLimMode', 'manual');
[xo3, yo3] = offsetCurve(x, y, 6, gca);
plot(xo3, yo3, '-m');

hold off;
