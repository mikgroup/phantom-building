load /Users/jtamir/Desktop/2016-09-18-phantom_test/T2fit_20160918.mat
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
slices = [5];
% idx2 = [2, 5, 9, 1, 6, 10, 3, 7, 11, 4, 8, 12];
% idx3 = [2, 6, 9, 1, 5, 10, 3, 7, 11, 4, 8, 12];
% idx4 = [2, 6, 9, 1, 5, 10, 3, 7, 11, 4, 8, 12];
% idx5 = [2, 5, 9, 1, 6, 10, 3, 7, 11, 4, 8, 12];
% idx6 = [2, 5, 9, 1, 6, 10, 3, 7, 11, 4, 8, 12];
% idx10 = [2, 5, 9, 1, 6, 10, 3, 7, 11, 4, 8, 12];
% idx11 = [2, 5, 9, 1, 6, 10, 3, 7, 11, 4, 8, 12];
% idx12 = [2, 6, 9, 1, 5, 10, 3, 7, 11, 4, 8, 12];
% idx13 = [2, 6, 9, 1, 5, 10, 3, 7, 11, 4, 8, 12];
% idx14 = [2, 5, 9, 1, 6, 10, 3, 7, 11, 4, 8, 12];
idx5 = [3, 6, 2, 5, 1, 4];

idxs = {idx5};

R1vals = cell(1, length(slices));

map = T2est;
VV = {};
for ii=1:length(slices)
    sl = slices(ii);
    idx = idxs{ii};
    x1 = squeeze(map(:,:,sl,:));
    
    v = zeros(length(idx), 1);
    v2 = cell(length(idx), 1);
    for jj=1:2:length(idxs{ii})
        m1 = (labels(:,:,sl)==idx(jj)) + (labels(:,:,sl)==idx(jj+1));
        x2 = sort(1 ./ x1(repmat(m1, [1, 1, size(x1,3)])==1));
        i1 = find(x2 > .01*median(x2), 1, 'first');
        i2 = find(x2 <= .99*median(x2), 1, 'last');
        x3 = x2(i1:i2); % throw out outliers
        v(jj) = mean(x3);
        v2{jj} = x3;
    end
    v = v(1:2:end);
    v2 = v2(1:2:end);
    VV{end+1} = v2;
    
    R1vals{ii} = reverse(v);
end

%%
axis1 = [.3, .6, .9]; % mM MnCl2

xlabel1 = 'mM MnCl2';

xlabels = {xlabel1};

axis = {axis1};
slices = 1;
R1trend = zeros(2, length(slices));

for ii=1:length(slices)
    [P, S] = polyfit(axis{ii}.', R1vals{ii}, 1);
    [Y, E] = polyconf(P, axis{ii}, S);
    R1trend(:, ii) = P;
    
    figure(ii);
    plot(axis{ii}, R1vals{ii}, 'o', 'linewidth', 3)
    hold on;
    errorbar(axis{ii}, Y, E, 'k--', 'linewidth', 2);
    hold off
    xlabel(xlabels{ii});
    ylabel('R1 (1/s)');
    legend('data', 'fit with 95% CI');
    
    fprintf('%20s:\tm1 = %.3f\tR1w = %.3f\n', xlabels{ii}, R1trend(1, ii), R1trend(2, ii));
end
