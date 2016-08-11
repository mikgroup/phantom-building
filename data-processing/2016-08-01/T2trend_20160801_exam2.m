load ~/Desktop/2016-08-01_phantom-test/exam2/T2fit.mat
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
slices = [30, 37, 37, 37];
idx30 = [1, 6, 10, 2, 5, 9, 3, 7, 11, 4, 8, 12];
idx37 = [5, 14, 21, 28, 36, 2, 10, 17, 23, 30, 1, 13];
idx37_2 = [20, 27, 34, 4, 12, 16, 26, 33, 3, 9, 19, 25]; 
idx37_3 = [32, 7, 11, 18, 24, 31, 6, 15, 22, 29, 35, 8];
idxs = {idx30, idx37, idx37_2, idx37_3};

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
    
    R2vals{ii} = v;
end

%%
axis30 = [6, 8, 10, 12, 14, 16]; % mM NiCl2
axis37 = [6, 8, 10, 12, 14, 16]/5; % mM NiCl2
axis37_2 = [6, 8, 10, 12, 14, 16]/5; % mM Co(NO3)2
axis37_3 = [6, 8, 10, 12, 14, 16]/5; % mM CuSO4

xlabel30 = 'mM NiCl2 (good)';
xlabel37 = 'mM NiCl2 (bad)';
xlabel37_2 = 'mM Co(NO3)2 (bad)';
xlabel37_3 = 'mM CuSO4 (bad)';

xlabels = {xlabel30, xlabel37, xlabel37_2, xlabel37_3};

axis = {axis30, axis37, axis37_2, axis37_3};

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
