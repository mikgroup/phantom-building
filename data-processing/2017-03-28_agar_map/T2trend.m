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
idx1 = flip([7, 8, 10, 9, 12, 11, 4, 6, 3, 5, 2, 1]);

idxs = {idx1};

R2vals = cell(1, length(slices));

map = T2est;

for sl = 1:length(slices)
    labels(:,:,sl) = labels(:,:,sl) .* flip(flip(m,2),1);
end

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
        i2 = find(x2 <= 1.99*median(x2), 1, 'last');
        x3 = x2(i1:i2); % throw out outliers
        v(jj) = mean(x3);
        v2{jj} = x3;
    end
    
    R2vals{ii} = v(end:-1:1);
    v2 = flip(v2);
end

%%

targetT1 = [830, 500, 1230, 1000, 250, 500];
targetT2 = [80, 110, 110, 120, 40, 60];
baseWV = 4.5 / 300;
vialTotal = 20;
vialGelVol = [0, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 20];
vialWV = vialGelVol/vialTotal * baseWV;

numPoints = 0;
for i = 1:length(v2)
    numPoints = numPoints + length(v2{i});
end

yMat = zeros(numPoints,1);
A = zeros(numPoints, 4);
ind = 1;


plot(vialWV, 1000./R2vals{1});

%A \ yMat