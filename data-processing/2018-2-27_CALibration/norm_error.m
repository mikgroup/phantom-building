close all
figure
plot(abs((1000./r1mix - t1_vals')./(1000./r1mix))./...
    max(abs((1000./r1mix - t1_vals')./(1000./r1mix))), '-o','linewidth',3,'MarkerSize',6)
hold on
%plot(abs((1000./r2mix - t2_vals')./(1000./r2mix))./...
%    max(abs((1000./r2mix - t2_vals')./(1000./r2mix))), '-o','linewidth',3,'MarkerSize',6)
plot(kAmix./max(kAmix),'-o','linewidth',3,'MarkerSize',6)
%plot(kBmix./max(kBmix),'-o','linewidth',3,'MarkerSize',6)
grid on
legend('T1 Error', 'Volume')
xlabel('Sample #', 'FontSize', 20)
ylabel('Normalized Values', 'FontSize', 20)

