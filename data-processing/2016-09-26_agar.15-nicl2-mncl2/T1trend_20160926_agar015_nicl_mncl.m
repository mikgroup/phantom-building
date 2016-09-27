load T1fit.mat
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
idx1 = [2, 7, 9, 1, 6, 11]; % NiCl2
idx2 = [4, 5, 10, 3, 8, 12]; % MnCl2
idxs = {idx1, idx2};

R1vals = cell(1, length(slices));

map = T1est;

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
    
    R1vals{ii} = v(end:-1:1);
end

%%
axis1 = [6, 8, 10, 12, 14, 16]; % mM NiCl2
axis2 = [.05, .1, .15, .2, .25, .3]; % mM MnCl2


xlabel1 = 'mM NiCl2';
xlabel2 = 'mM MnCl2';

xlabels = {xlabel1, xlabel2};

axiss = {axis1, axis2};

R1trend = zeros(2, length(slices));

for ii=1:length(slices)
    [P, S] = polyfit(axiss{ii}.', R1vals{ii}, 1);
    [Y, E] = polyconf(P, axiss{ii}, S);
    R1trend(:, ii) = P;
    
    figure(ii);
    plot(axiss{ii}, R1vals{ii}, 'o', 'linewidth', 3)
    hold on;
    errorbar(axiss{ii}, Y, E, 'k--', 'linewidth', 2);
    hold off
    xlabel(xlabels{ii});
    ylabel('R1 (1/s)');
    legend('data', 'fit with 95% CI');
    axis square;
    fprintf('%20s:\tm1 = %.3f\tR1w = %.3f\n', xlabels{ii}, R1trend(1, ii), R1trend(2, ii));
end
