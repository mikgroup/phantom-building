bigA = [];
bigY = [];

T1trend;
At1 = A;
yt1 = yMat;
t1Vals = 1000./v;

T2trend;
At2 = A;
yt2 = yMat;
bigA = [At1 ; At2];
bigY = [yt1 ; yt2];
t2Vals = 1000./v;

clearvars -except bigA bigY t1Vals t2Vals
close all

sol = bigA \ bigY;

% niclConc = [0,1.5775,0,0.8048,0,0.7255];
% mnclConc = [0.4449,0.3146,0.7077,0.6412,0.5494,0.4895];
% agarConc = [0.0054,0.0075,0.0064,0.0075,0.0165,0.0175];
% 
% res = zeros(length(niclConc), 2);
% for i = 1:length(niclConc)
%     res(i,:) = getT1T2FromConc(sol, niclConc(i), mnclConc(i), agarConc(i));
% end

%%
%% Paramagnetic Ions
% Ni  Mn
%
ions = {'NiCl_2', 'MnCl_2', 'Agar'};
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
legend({sprintf('%s bound', ions{ion1}), sprintf('%s bound', ions{ion2}), sprintf('%s bound', ions{ion3}), 'target values'}, 'Location', 'northwest')
faxis(gca, 20);
hold off

figure(2);

t2Tissue = [42,50,47,56,27,43,69,99,78,78,275];
t1Tissue = [812,1412,1471,1194,1168,1156,1084,1820,1083,993,1932];

hold on
plot(1000./r2_range, 1000./r2fun(ion1),'-', 1000./r2_range, 1000./r2fun(ion2),'-', 1000./r2_range, 1000./r2fun(ion3),'-', t2Vals, t1Vals, '+', 'linewidth', 3, 'MarkerSize', 10)

grid on
grid minor
xlabel('T_2 (ms)');
ylabel('T_1 (ms)');
xlim([0, 500]);
ylim([0, 5000]);
legend({sprintf('%s bound', ions{ion1}), sprintf('%s bound', ions{ion2}),sprintf('%s bound', ions{ion3}), 'Data'}, 'Location', 'northwest');
faxis(gca, 20);
set(0,'DefaultAxesLineStyleOrder',{'+','o','*','.','x','s','d','^','v','>','<','p','h'});
for i = 1:length(t2Tissue)
    plot(t2Tissue(i),t1Tissue(i), 'MarkerSize', 10)
end
%plot(t2_vals,t1_vals,'ro');
saveas(2, '~/Google Drive Berkeley/phantom-building/figs/10_19_mapping/T1_T2_Bounds.svg');

hold off

% Concentrations in mM
fprintf('mM Concentrations of %s\n',ions{ion1})
kAmix = (m2(ion2) * (r1mix - rW1fixed) - m1(ion2) * (r2mix - rW2fixed))./(m2(ion2) * m1(ion1) - m1(ion2) * m2(ion1))
fprintf('Max concentration of %s: %.03f mM\n', ions{ion1}, max(kAmix));
kBmix = -(m2(ion1) * (r1mix - rW1fixed) - m1(ion1) * (r2mix - rW2fixed))./(m2(ion2) * m1(ion1) - m1(ion2) * m2(ion1))
fprintf('Max concentration of %s: %.03f mM\n', ions{ion2}, max(kBmix));