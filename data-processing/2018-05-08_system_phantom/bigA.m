bigAmat = [];
bigYmat = [];

T1trend;
bigAmat = [bigAmat ; A];
bigYmat = [bigYmat ; yMat];
t1_vals = 1000./v;
clearvars -except bigAmat bigYmat t1_vals

T2trend;
bigAmat = [bigAmat ; A];
bigYmat = [bigYmat ; yMat];
t2_vals = 1000./v;
clearvars -except bigAmat bigYmat t1_vals t2_vals

close all

sol = bigAmat \ bigYmat;

% t1_vals = [756, 522, 373,283,380,1102,823,528,923,1026,833,1446,1261,1076,1314,1498];
% t2_vals = [70.0, 54.4, 57.5, 78.0, 97.7, 94.9, 99.5, 108.4, 109.7, 116.6,...
% 140.5, 120.4, 133.2, 144.8, 151.7, 155.3];
% 

%% Paramagnetic Ions
% Ni  Mn
%
ions = {'NiCl_2', 'MnCl_2'};
molarMass = [129.6, 125.844];
m1 = [sol(1), sol(2)];
rW1 = [sol(5), sol(5)];
m2 = [sol(3), sol(4)];
rW2 = [sol(6), sol(6)];

ion1 = 1;
ion2 = 2;
fprintf('Ion 1 = %s, Ion 2 = %s\n', ions{ion1}, ions{ion2})

r2_range = 0.001:0.01:1000;
rW1fixed = mean(rW1);
rW2fixed = mean(rW2);

r2fun = @(ii) m1(ii) / m2(ii) * (r2_range - rW2fixed) + rW1fixed;
t2fun = @(ii) m2(ii) / (m1(ii)*(r2_range - rW2fixed) + rW1fixed);


%t2mix = [80, 110, 110, 120, 40, 60];
%t1mix = [830, 500, 1230, 1000, 250, 500];

%t2mix = [99, 70, 130];
%t1mix = [1290, 930, 1500];

%t2mix = [120.621545744990,120.644811528621,114.072523786058,119.407695074628];
%t1mix = [1700.19986775720,1718.92090125473,1630.71600808535,1621.97490069183];

% Cal letters
t1_C = [800,549.052388289677,387.912172573190,288.158705701079,387.912172573190];
t1_A = [1116.87981510015,848.312788906010,541.379044684129,932.719568567027,1047.81972265023,863.659476117103];
t1_L = [1508.22033898305,1293.36671802774,1109.20647149461,1347.08012326656,1546.58705701079];

t2_C = [65,52.9488703923900,56.9441141498217,77.4910820451843,96.8965517241380];
t2_A = [92.3305588585018,97.4673008323425,108.311533888228,110.594530321046,120.297265160523,145.980975029727];
t2_L = [120.297265160523,136.278240190250,153.400713436385,164.244946492271,172.235434007134];

t1mix = [t1_C,t1_A,t1_L];
t2mix = [t2_C,t2_A,t2_L];


r2mix = 1000./t2mix;
r1mix = 1000./t1mix;

figure(1);
plot(r2_range, r2fun(ion1), r2_range, r2fun(ion2), r2_range, 0 * r2_range, 'k--', r2mix, r1mix, '+', 'linewidth', 3, 'MarkerSize', 10)
grid on
grid minor
xlabel('R_2 (1/s)');
ylabel('R_1 (1/s)');
xlim([4, 30]);
ylim([.1, 15]);
legend({sprintf('%s bound', ions{ion1}), sprintf('%s bound', ions{ion2}), 'target values'}, 'Location', 'northwest')
faxis(gca, 20);

figure(2);
plot(1000./r2_range, 1000./r2fun(ion1), 1000./r2_range, 1000./r2fun(ion2), 1000./r2mix, 1000./r1mix, '+', 'linewidth', 3, 'MarkerSize', 10)
hold on
plot(t2_vals, t1_vals, '+', 'linewidth', 3, 'MarkerSize', 10)
hold off
grid on
grid minor
xlabel('T_2 (ms)');
ylabel('T_1 (ms)');
xlim([10, 250]);
ylim([10, 2500]);
legend({sprintf('%s bound', ions{ion1}), sprintf('%s bound', ions{ion2}), 'Target Values'}, 'Location', 'northwest');
faxis(gca, 20);
saveas(2, '~/Google Drive Berkeley/phantom-building/figs/10_19_mapping/T1_T2_Bounds.svg');


% Concentrations in mM
fprintf('mM Concentrations of %s\n',ions{ion1})
kAmix = (m2(ion2) * (r1mix - rW1fixed) - m1(ion2) * (r2mix - rW2fixed))./(m2(ion2) * m1(ion1) - m1(ion2) * m2(ion1))
fprintf('Max concentration of %s: %.03f mM\n', ions{ion1}, max(kAmix));
kBmix = -(m2(ion1) * (r1mix - rW1fixed) - m1(ion1) * (r2mix - rW2fixed))./(m2(ion2) * m1(ion1) - m1(ion2) * m2(ion1))
fprintf('Max concentration of %s: %.03f mM\n', ions{ion2}, max(kBmix));

hold off