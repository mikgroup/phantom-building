load ~/Desktop/2016-09-07_phantom-test/T2fit.mat
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
% slices = [17, 23];
% idx17 = [4, 7, 11, 3, 8, 12, 5, 9, 13, 6, 10, 14];
% idx23 = [7, 7, 12, 4, 6, 11, 3, 9, 13, 5, 8, 10]; % repeat vial 7 due to aliasing error
% idxs = {idx17, idx23};

slices = [1];
idx1 = [2, 6, 9, 1, 5, 11, 4, 8, 10, 3, 7, 12];
idxs = {idx1};

R2vals = cell(1, length(slices));

map = T2est;

for ii=1:length(slices)
    sl = slices(ii);
    idx = idxs{ii};
    x1 = squeeze(map(:,:,sl,:));
    
    v = zeros(length(idx), 1);
    v2 = cell(length(idx), 1);
    for jj=1:2:length(idxs{ii})
        m1 = (labels(:,:,sl)==idx(jj)) + (labels(:,:,sl)==idx(jj+1)) ~= 0;
        x2 = sort(1 ./ x1(repmat(m1, [1, 1, size(x1,3)])==1));
        i1 = find(x2 > .01*median(x2), 1, 'first');
        i2 = find(x2 <= .99*median(x2), 1, 'last');
        x3 = x2(i1:i2); % throw out outliers
        v(jj) = mean(x3);
        v2{jj} = x3;
    end
    v = v(1:2:end);
    v2 = v2(1:2:end);
    
    R2vals{ii} = reverse(v);
end

%%
axis1 = [.2, .4, .6, .8, 1, 1.2]; % mM MnCl2


xlabel1 = 'mM MnCl2';

xlabels = {xlabel1};

axis = {axis1};

R2trend = zeros(2, length(slices));

for ii=1:length(slices)
    [P, S] = polyfit(axis{ii}.', R2vals{ii}, 1);
    [Y, E] = polyconf(P, axis{ii}, S);
    R2trend(:, ii) = P;
    
    figure(ii);
    plot(axis{ii}, R2vals{ii}, 'o', 'linewidth', 3)
    hold on;
    errorbar(axis{ii}, Y, E, 'k--', 'linewidth', 2);
    hold off
    xlabel(xlabels{ii});
    ylabel('R2 (1/s)');
    legend('data', 'fit with 95% CI');
    
    fprintf('%20s:\tm2 = %.3f\tR2w = %.3f\n', xlabels{ii}, R2trend(1, ii), R2trend(2, ii));
end
