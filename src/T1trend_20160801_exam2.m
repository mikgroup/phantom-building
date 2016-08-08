load ~/Desktop/2016-08-01_phantom-test/exam2/T1fit.mat
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
slices = [32, 37, 37, 37];
idx32 = [2, 7, 10, 4, 6, 11, 3, 8, 13, 5, 9, 12];
idx37 = [8, 15, 22, 30, 37, 2, 10, 17, 24, 31, 5, 14];
idx37_2 = [20, 28, 35, 4, 13, 19, 27, 34, 3, 12, 18, 26];
idx37_3 = [33, 7, 11, 21, 25, 32, 6, 16, 23, 29, 36, 9];
idxs = {idx32, idx37, idx37_2, idx37_3};

R1vals = cell(1, length(slices));

map = T1est;

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
    
    R1vals{ii} = v;
end

%%
axis32 = [6, 8, 10, 12, 14, 16]; % mM NiCl2
axis37 = [6, 8, 10, 12, 14, 16]/5; % mM NiCl2
axis37_2 = [6, 8, 10, 12, 14, 16]/5; % mM Co(NO3)2
axis37_3 = [6, 8, 10, 12, 14, 16]/5; % mM CuSO4

xlabel32 = 'mM NiCl2 (good)';
xlabel37 = 'mM NiCl2 (bad)';
xlabel37_2 = 'mM Co(NO3)2 (bad)';
xlabel37_3 = 'mM CuSO4 (bad)';

xlabels = {xlabel32, xlabel37, xlabel37_2, xlabel37_3};

axis = {axis32, axis37, axis37_2, axis37_3};

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
