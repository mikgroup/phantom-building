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
t1_C = [809.946070878274,549.052388289677,387.912172573190,288.158705701079,387.912172573190];
t1_A = [1116.87981510015,848.312788906010,541.379044684129,932.719568567027,1047.81972265023,863.659476117103];
t1_L = [1508.22033898305,1293.36671802774,1109.20647149461,1347.08012326656,1546.58705701079];

t2_C = [62.6516052318668,52.9488703923900,56.9441141498217,77.4910820451843,96.8965517241380];
t2_A = [92.3305588585018,97.4673008323425,108.311533888228,110.594530321046,120.297265160523,145.980975029727];
t2_L = [120.297265160523,136.278240190250,153.400713436385,164.244946492271,172.235434007134];

%t1mix = [t1_C,t1_A,t1_L];
%t2mix = [t2_C,t2_A,t2_L];

t1mix = [111.671802773498,272.812018489985,257.465331278891,142.365177195686,142.365177195686,357.218798151002,495.338983050848,495.338983050848,334.198767334361,295.832049306626,526.032357473036,641.132511556240,763.906009244992,909.699537750385,863.659476117103,702.519260400617,449.298921417566,556.725731895224,855.986132511556,925.046224961479,1116.87981510015,1331.73343605547,1431.48690292758,1377.77349768875,1201.28659476117,1139.89984591680,1324.06009244992,1416.14021571649,1584.95377503852,1799.80739599384,1983.96764252696,1815.15408320493,1822.82742681048,2091.39445300462,2183.47457627119];
t2mix = [30.1189060642093,24.4114149821641,38.1093935790725,44.9583828775268,60.9393579072533,53.5196195005945,42.6753864447087,72.9250891795482,76.3495838287753,90.0475624256837,91.1890606420928,84.3400713436385,62.6516052318668,76.3495838287753,90.0475624256837,111.736028537455,120.297265160523,133.424494649227,131.141498216409,111.165279429251,106.028537455410,100.891795481570,115.731272294887,130.570749108205,141.985731272295,163.103448275862,174.518430439952,160.820451843044,150.546967895363,138.561236623068,152.259215219976,165.957193816885,183.650416171225,183.650416171225,166.527942925089];

samp_sol_1 = [1.8279,0.2359,2.7895,1.1975,0.4127,3.2427,1.6507,0.8659,0.3985,0.0883,1.9144,1.1295,0.6621,0.352,0.1312,1.302,0.8346,0.5245,0.3036,0.1384,0.9562,0.6461,0.4253,0.26,0.1317,0.3504,0.2221,0.1196];
samp_sol_2 = [0.5625,0.6246,0.2381,0.3002,0.3308,0.0852,0.1473,0.1779,0.1961,0.2082,0.0583,0.0889,0.1071,0.1192,0.1279,0.0307,0.049,0.0611,0.0697,0.0761,0.0079,0.02,0.0286,0.0351,0.0401,0.0046,0.0096,0.0136];


t1mix = (1000./(m1(1) * samp_sol_1*1.25 + m1(2) .* samp_sol_2*.2 + rW1fixed));
t2mix = (1000./(m2(1) * samp_sol_1*1.25 + m2(2) .* samp_sol_2*.2 + rW2fixed));


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
plot(1000./r2_range, 1000./r2fun(ion1), 1000./r2_range, 1000./r2fun(ion2), 'linewidth', 3, 'MarkerSize', 10)
hold on
plot(t2_vals, t1_vals, '+', 'linewidth', 3, 'MarkerSize', 10)
hold off
grid on
grid minor
xlabel('T_2 (ms)');
ylabel('T_1 (ms)');
xlim([10, 250]);
ylim([10, 2500]);
legend({sprintf('%s bound', ions{ion1}), sprintf('%s bound', ions{ion2}), 'Data'}, 'Location', 'northwest');
faxis(gca, 20);
%saveas(2, '~/Google Drive Berkeley/phantom-building/figs/10_19_mapping/T1_T2_Bounds.svg');

figure(3);
plot(1000./r2_range, 1000./r2fun(ion1), 1000./r2_range, 1000./r2fun(ion2), 'linewidth', 3, 'MarkerSize', 10)
hold on
plot(t2mix, t1mix, '+', 'linewidth', 3, 'MarkerSize', 10)
hold off
grid on
grid minor
xlabel('T_2 (ms)');
ylabel('T_1 (ms)');
xlim([10, 250]);
ylim([10, 2500]);
legend({sprintf('%s bound', ions{ion1}), sprintf('%s bound', ions{ion2}), 'Sampling'}, 'Location', 'northwest');
faxis(gca, 20);
%saveas(2, '~/Google Drive Berkeley/phantom-building/figs/10_19_mapping/T1_T2_Bounds.svg');


% Concentrations in mM
fprintf('mM Concentrations of %s\n',ions{ion1})
kAmix = (m2(ion2) * (r1mix - rW1fixed) - m1(ion2) * (r2mix - rW2fixed))./(m2(ion2) * m1(ion1) - m1(ion2) * m2(ion1))
fprintf('Max concentration of %s: %.03f mM\n', ions{ion1}, max(kAmix));
kBmix = -(m2(ion1) * (r1mix - rW1fixed) - m1(ion1) * (r2mix - rW2fixed))./(m2(ion2) * m1(ion1) - m1(ion2) * m2(ion1))
fprintf('Max concentration of %s: %.03f mM\n', ions{ion2}, max(kBmix));

hold off