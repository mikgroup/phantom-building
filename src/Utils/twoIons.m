%% Paramagnetic Ions
% Ni Cu Co Mn
clear;
close all;
ions = {'NiCl_2', 'CuSO_4', 'Co(NO_3)_2', 'MnCl_2'};
molarMass = [129.6 159.61 291.03 125.844];
m1 = [0.567 0.746 0.122 9.804];
rW1 = [0.703 1.161 1.098 -0.301];
m2 = [1.230 1.669 0.343 76.577];
rW2 = [7.241 11.015 11.669 12.105];

div = m1 ./ m2;
divMat = [div ; div ; div ; div];

subRatios = divMat - divMat';

maxDiff = max(subRatios(:));
[ion1, ion2] = find(subRatios == maxDiff);
%ion1 = 2;
%ion2 = 3;
fprintf('Ion 1 = %s, Ion 2 = %s\n', ions{ion1}, ions{ion2})

t = 0:0.01:30;
rW1fixed = (1.161 + 1.098)/2;
rW2fixed = (11.015 + 11.669)/2;

figure
xlabel('R_2');
ylabel('R_1');
hold on
plot(t, div(ion1) * (t - rW2fixed) + rW1fixed)
plot(t, div(ion2) * (t - rW2fixed) + rW1fixed)
plot(t, 0 * t)

plot([25 20 12 30 15 25], [3 3 1.25 9.4 2.5 6], 'o')
%plot([30 30 26 15 20 12], [8 9.2 7.5 2.5 4.8 1.4], 'o')
grid on
grid minor

%%
R2mix1 = [25 20 12 30 15 25];
R1mix1 = [3 3 1.25 9.4 2.5 6];

% Concentrations in mM
fprintf('mM Concentrations of %s',ions{ion1})
kAmix1 = (m2(ion2) * (R1mix1 - rW1fixed) - m1(ion2) * (R2mix1 - rW2fixed))./(m2(ion2) * m1(ion1) - m1(ion2) * m2(ion1))
max(kAmix1)
fprintf('mM Concentrations of %s',ions{ion2})
kBmix1 = -(m2(ion1) * (R1mix1 - rW1fixed) - m1(ion1) * (R2mix1 - rW2fixed))./(m2(ion2) * m1(ion1) - m1(ion2) * m2(ion1))
max(kBmix1)

% ion1 = 2;
% ion2 = 3;
% R2mix2 = [30 30 26 15 20 12];
% R1mix2 = [8 9.2 7.5 2.5 4.8 1.4];
% 
% kAmix2 = (m2(ion2) * (R1mix2 - rW1fixed) - m1(ion2) * (R2mix2 - rW2fixed))./(m2(ion2) * m1(ion1) - m1(ion2) * m2(ion1))
% kBmix2 = -(m2(ion1) * (R1mix2 - rW1fixed) - m1(ion1) * (R2mix2 - rW2fixed))./(m2(ion2) * m1(ion1) - m1(ion2) * m2(ion1))


hold off
