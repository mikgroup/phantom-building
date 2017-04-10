%%
% Produces graph for 3 values of Agar
close all;
clear;
bigAmat = [];
bigYmat = [];

t1_vals = {};
t2_vals = {};
sol = {};

T1trend_low_agar;
bigAmat = [bigAmat ; A];
bigYmat = [bigYmat ; yMat];
t1_vals{1} = 1000./v;
clearvars -except bigAmat bigYmat

T2trend_low_agar;
bigAmat = [bigAmat ; A];
bigYmat = [bigYmat ; yMat];
t2_vals{1} = 1000./v;
sol{1} = bigAmat \ bigYmat;
clearvars -except bigAmat bigYmat


T1trend_medium_agar;
bigAmat = [bigAmat ; A];
bigYmat = [bigYmat ; yMat];
t1_vals{2} = 1000./v;
clearvars -except  bigAmat bigYmat

T2trend_medium_agar;
bigAmat = [bigAmat ; A];
bigYmat = [bigYmat ; yMat];
t2_vals{2} = 1000./v;
sol{2} = bigAmat \ bigYmat;
clearvars -except bigAmat bigYmat


T1trend_high_agar;
bigAmat = [bigAmat ; A];
bigYmat = [bigYmat ; yMat];
t1_vals{3} = 1000./v;
clearvars -except bigAmat bigYmat

T2trend_high_agar;
bigAmat = [bigAmat ; A];
bigYmat = [bigYmat ; yMat];
t2_vals{3} = 1000./v;
sol{3} = bigAmat \ bigYmat;
clearvars -except bigAmat bigYmat

sol = bigAmat \ bigYmat;
close all;


%% Paramagnetic Ions
% Ni  Mn
%
ions = {'NiCl_2', 'MnCl_2'};
molarMass = [129.6, 125.844];
m1 = [sol(1), sol(2), sol(3)];
rW1 = [sol(7), sol(7), sol(7)];
m2 = [sol(4), sol(5), sol(6)];
rW2 = [sol(8), sol(8), sol(8)];

ion1 = 1;
ion2 = 2;
ion3 = 3;
fprintf('Ion 1 = %s, Ion 2 = %s\n', ions{ion1}, ions{ion2})

r2_range = 0.001:0.01:1000;
rW1fixed = mean(rW1);
rW2fixed = mean(rW2);

r2fun = @(ii) m1(ii) / m2(ii) * (r2_range - rW2fixed) + rW1fixed;
t2fun = @(ii) m2(ii) / (m1(ii)*(r2_range - rW2fixed) + rW1fixed);


t2mix = [137.0454  155.7443  151.3693  156.5790   90.7592  116.6885];
t1mix = [1.9183    1.3135    2.3473    2.1285    0.7584    1.3140] * 1e3;

% t2mix = [99, 70];
% t1mix = [1290, 930];

r2mix = 1000./t2mix;
r1mix = 1000./t1mix;

figure(1);
hold on
plot(r2_range, r2fun(ion1), r2_range, r2fun(ion2), r2_range, r2fun(ion3), r2_range, 0 * r2_range, 'k--', r2mix, r1mix, '+', 'linewidth', 3, 'MarkerSize', 10)
grid on
grid minor
xlabel('R_2 (1/s)');
ylabel('R_1 (1/s)');
xlim([4, 30]);
ylim([.1, 15]);
legend({sprintf('%s bound', ions{ion1}), sprintf('%s bound', ions{ion2}), 'target values'}, 'Location', 'northwest')
faxis(gca, 20);
hold off

figure(2);
hold on
plot(1000./r2_range, 1000./r2fun(ion1), 1000./r2_range, 1000./r2fun(ion2), 1000./r2_range, 1000./r2fun(ion3), 1000./r2mix, 1000./r1mix, '+', 'linewidth', 3, 'MarkerSize', 10)
grid on
grid minor
xlabel('T_2 (ms)');
ylabel('T_1 (ms)');
xlim([10, 500]);
ylim([10, 5000]);
%legend({sprintf('%s bound', ions{ion1}), sprintf('%s bound', ions{ion2}), 'Target Values'}, 'Location', 'northwest');
faxis(gca, 20);
plot(t2_vals,t1_vals,'ro');
saveas(2, '~/Google Drive Berkeley/phantom-building/figs/10_19_mapping/T1_T2_Bounds.svg');

hold off

% Concentrations in mM
fprintf('mM Concentrations of %s\n',ions{ion1})
kAmix = (m2(ion2) * (r1mix - rW1fixed) - m1(ion2) * (r2mix - rW2fixed))./(m2(ion2) * m1(ion1) - m1(ion2) * m2(ion1))
fprintf('Max concentration of %s: %.03f mM\n', ions{ion1}, max(kAmix));
kBmix = -(m2(ion1) * (r1mix - rW1fixed) - m1(ion1) * (r2mix - rW2fixed))./(m2(ion2) * m1(ion1) - m1(ion2) * m2(ion1))
fprintf('Max concentration of %s: %.03f mM\n', ions{ion2}, max(kBmix));

