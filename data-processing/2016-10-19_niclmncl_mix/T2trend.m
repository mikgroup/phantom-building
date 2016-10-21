load T2fit.mat
rng(10);

mask = T2mask;

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
slices = [1];
idx1 = [2, 7, 5, 12, 1, 8, 4, 9, 3, 11, 6, 10];

idxs = {idx1};

R1vals = cell(1, length(slices));

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
    
    R1vals{ii} = v(end:-1:1);
    v2 = flip(v2);
end

%%

targetT1 = [830, 500, 1230, 1000, 250, 500];
targetT2 = [80, 110, 110, 120, 40, 60];
concA = [0.7725, 5.2487, 0.3245, 1.3619, 8.0507, 2.6321];
concB = [0.1715, 0.0453, 0.0986, 0.0701, 0.3776, 0.2463];

order = [3, 3, 4, 4, 1, 1, 2, 2, 6, 6, 5, 5];

axis1 = zeros(1, length(order));
axis2 = zeros(1, length(order));
for i = 1:length(order)
   axis1(i) = concA(order(i));
   axis2(i) = concB(order(i));
end

xlabel1 = 'mM NiCl2';
xlabel2 = 'mM MnCl2';

xlabels = {xlabel1, xlabel2};

axiss = {axis1, axis2};

numPoints = 0;
for i = 1:length(v2)
    numPoints = numPoints + length(v2{i});
end

yMat = zeros(numPoints,1);
A = zeros(numPoints, 6);
ind = 1;

t1 = false;

for i = 1:length(v2)
    kA = axis1(i);
    kB = axis2(i);
    
    for j = 1:length(v2{i})
       if (t1)
        A(ind,:) = [kA kB 0 0 1 0];
        yMat(ind) = v2{i}(j);
       else
        A(ind,:) = [0 0 kA kB 0 1];
        yMat(ind) = v2{i}(j);
       end
       ind = ind + 1;
    end
end


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

%A \ yMat
