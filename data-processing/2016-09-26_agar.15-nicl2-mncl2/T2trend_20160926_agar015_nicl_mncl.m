load T2fit.mat
rng(10);
labels_cc = zeros(size(mask));
SE = strel('diamond',2);
for ii=1:ns
    m0 = mask(:,:,ii);
    m1 = imerode(m0, SE);
    [L, ~] = bwlabel(m1, 8);
    labels_cc(:,:,ii) = L;
end

labels = labels_cc;
clear labels_cc m0 m1 L SE

num = max(reshape(labels, [], ns), [], 1).';

%%
slices = [5, 5];
idx1 = [4, 8, 9, 1, 7, 12]; % NiCl2
idx2 = [3, 6, 11, 2, 5, 10]; %MnCl2
idxs = {idx1, idx2};

R2vals = cell(1, length(slices));

map = T2est;

for ii=1:length(slices)
    sl = slices(ii);
    idx = idxs{ii};
    x1 = squeeze(map(:,:,sl,:));
    
    v = zeros(length(idx), 1);
    v2 = cell(length(idx), 1);
    for jj=1:length(idxs{ii})
        m1 = (labels(:,:,sl)==idx(jj));
        x2 = sort(1 ./ x1(repmat(m1, [1, 1, size(x1,3)])==1));
        i1 = find(x2 > .01*median(x2), 1, 'first');
        i2 = find(x2 <= .99*median(x2), 1, 'last');
        x3 = x2(i1:i2); % throw out outliers
        v(jj) = mean(x3);
        v2{jj} = x3;
    end
    
    R2vals{ii} = v(end:-1:1);
end

%%
axis1 = [6, 8, 10, 12, 14, 16]; % mM NiCl2
axis2 = [.05, .1, .15, .2, .25, .3]; % mM MnCl2


xlabel1 = 'NiCl_2 Concentration (mM)';
xlabel2 = 'MnCl_2 Concentration (mM)';

xlabels = {xlabel1, xlabel2};

figNames = {'NiCl2_T2', 'MnCl2_T2'};

axiss = {axis1, axis2};

R2trend = zeros(2, length(slices));

for ii=1:length(slices)
    [P, S] = polyfit(axiss{ii}.', R2vals{ii}, 1);
    [Y, E] = polyconf(P, axiss{ii}, S);
    extendPerc = 0.15;
    lowerExtend = min(axiss{ii}) - extendPerc * mean([min(axiss{ii}), max(axiss{ii})]);
    upperExtend = max(axiss{ii}) + extendPerc * mean([min(axiss{ii}), max(axiss{ii})]);
    axisExtend = lowerExtend:0.01:upperExtend;
    [yExtend, eExtend] = polyconf(P, axisExtend, S);
    R2trend(:, ii) = P;
    
    fig(ii) = figure(ii*10);
    ax = gca;
    h = plot(axisExtend, yExtend, 'linewidth', 2);
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on;
    ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
    errorbar(axiss{ii}, Y, E, 'linewidth', 2, 'Capsize', 10);
    plot(axiss{ii}, R2vals{ii}, '+', 'MarkerSize', 10, 'LineWidth', 3)
    hold off
    xlabel(xlabels{ii});
    ylabel('R2 (1/s)');
    legend({'Fit (95% CI)', 'Data'}, 'Location', 'northwest');
    %xlim([0,0.35]);
    %ylim([0, 35]);
    axis tight
    axis square
    fprintf('%20s:\tm2 = %.3f\tR2w = %.3f\n', xlabels{ii}, R2trend(1, ii), R2trend(2, ii));
    faxis(gca, 20)
    saveas(fig(ii), strcat('~/Google Drive Berkeley/phantom-building/figs/9_27_mapping/',figNames{ii}), 'svg');
end
