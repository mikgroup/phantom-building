load ~/Desktop/2016-08-01_phantom-test/exam1/T1fit.mat

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
slices = [4, 10];

idx4 = [3, 5, 9, 2, 6, 10, 1, 7, 11, 4, 8, 12];
idx10 = [1, 6, 11, 2, 5, 10, 3, 8, 12, 4, 7, 9]; % repeat vial 7 twice since there is segmentation error with vial 1
idxs = {idx4, idx10};

R1vals = cell(1, length(slices));

map = T1est;

for ii=1:length(slices)
    sl = slices(ii);
    v = zeros(num(sl), 1);
    x1 = mean(map(:,:,sl,:),4);
    for jj=1:num(sl)
        v(jj) = mean(x1(labels(:,:,sl)==jj));
    end
    v = v(idxs{ii});
    
    R1vals{ii} = fliplr(1./mean(reshape(v, 2, []), 1));
end

%%
axis4 = [5, 20, 35, 50, 65, 80]; % mM Co(NO3)2
axis10 = [1, 2.8, 4.6, 6.4, 8.2, 10]; % mM CuSO4

xlabel4 = 'mM Co(NO3)2';
xlabel10 = 'mM CuSO4';
xlabels = {xlabel4, xlabel10};

axis = {axis4, axis10}; % FIXME: check if this is flipped!

R1trend = zeros(2, length(slices));

for ii=1:length(slices)
    P = polyfit(axis{ii}, R1vals{ii}, 1);
    R1trend(:, ii) = P;
    
    figure(ii);
    plot(axis{ii}, R1vals{ii}, '-o', axis{ii}, axis{ii}*R1trend(1,ii) + R1trend(2,ii), 'k--', 'linewidth', 3)
    xlabel(xlabels{ii});
    ylabel('R1 (1/s)');
    legend('data', 'fit');
end
