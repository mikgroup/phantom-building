close all

load T2fit.mat
rng(10);

mask = T2mask;
mask(T2est > 0.7) = 0;
mask(T2est < 0.01) = 0;

load('T1fit.mat', 'mask');

mask = (flip(mask,1));
T2est = (flip(T2est,1));

labels_cc = zeros(size(mask));
boundaries = cell(length(ns), 1);
SE = strel('diamond',3);
for ii=1:ns
    m0 = mask(:,:,ii);
    m1 = double(imerode(m0, SE));
    m1BW = imbinarize(m1);
    m1BW = bwareaopen(m1BW, 200); % Remove small regions
    
    % Labels and creates cell array of boundaries
    [B, L] = bwboundaries(m1BW, 'noholes'); 
    labels_cc(:,:,ii) = L;
    boundaries{ii} = B;
end

labels = labels_cc;
clear labels_cc m0 m1 L B SE

%%

for sl = 1:size(labels,3)
    % Plots boundaries
    % make sure these are not touching
    imagesc(label2rgb(labels(:,:,sl), @jet, [.5 .5 .5]))
    
    bl = 6;
    
    hold on
    for i = 1:length(boundaries{sl})
        boundary = boundaries{sl}{i};
        plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2);
    end
    
    % Find centers of blobs
    stats = regionprops(labels(:,:,sl), 'centroid');
    centers = cat(1,stats.Centroid);
    plot(centers(:,1), centers(:,2), 'ro');
    
    % Compute minimum inscribing circle for each blob
    radii = zeros(length(centers(:,1)),1);
    for blob = 1:length(centers(:,1))
        boundary = boundaries{sl}{blob};
        dists = zeros(length(boundary),1);
        for i = 1:length(boundary)
            dists(i) = pdist([flip(boundary(i,:)); centers(blob,:)]);
        end
        % Reduce radii to avoid boundary
        radii(blob) = 0.95 * min(dists);
        %viscircles(centers(blob,:), radii(blob))
    end
    
    [x, y] = meshgrid(1:length(labels(1,:,sl)), 1:length(labels(:,1,sl)));
    % Replace with radii to use different inscribing circles
    radiiToUse = radii;%repmat(min(radii), [size(radii), 1]);
    circMask = zeros(size(labels(:,:,sl)));
    
    for blob = 1:length(centers(:,1))
        viscircles(centers(blob,:), radiiToUse(blob), 'Color', 'b');
        
        circMask = circMask + ...
            double(sqrt((x - centers(blob,1)).^2 + ...
            (y - centers(blob,2)).^2) <= radiiToUse(blob));
        
        row = ceil(centers(blob,1));
        col = ceil(centers(blob,2));
        h = text(row, col, num2str(labels(col, row, 1)));
        set(h,'FontSize',14,'FontWeight','bold');
    end
    
    labels(:,:,sl) = labels(:,:,sl) .* circMask;
    hold off 
end
num = max(reshape(labels, [], ns), [], 1).';

%%
slices = [1];

idx1 = [5,10,16,19,2,11,13,20,3,7,17,21,1,8,14,23,4,9,15,22,6,12,18,24];
idx1 = order_labels(centers);
idxs = {idx1};

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
        i2 = find(x2 <= 1.99*median(x2), 1, 'last');
        x3 = x2(i1:i2); % throw out outliers
        v(jj) = mean(x3);
        v2{jj} = x3;
    end
    
    R2vals{ii} = v;
end

%%

%%
arrangement = ([2, 1, 8, 8, 7, 6, 1, 2, 5, 5, 4, 4, 7, 3, 3, 6]);

% nicl_sol_volume = [0.1068780411, 2.256814843, 3.298934513, 2.099942248, 0.689791923, 1.164271975, 0.3805046478, 0.5031443664];
% mncl_sol_volume = [0.09166293959, 0.09159891125, 0.01751858376, 0.01868012418, 0.02976340982, 0.005897133059, 0.01735628546, 0.005091547006];

nicl_sol_volume = [0.0855, 1.8055, 2.6391, 1.6800, 0.5518, 0.9314, 0.3044, 0.4025];
mncl_sol_volume = [0.3667, 0.3664, 0.0701, 0.0747, 0.1191, 0.0236, 0.0694, 0.0204];
agar = 15 * ones(1,16);

nicl_sol_volume = nicl_sol_volume(arrangement);
mncl_sol_volume = mncl_sol_volume(arrangement);
%agar = agar(arrangement);

nicl_conc = nicl_sol_volume ./ (agar + 5);

mncl_conc = mncl_sol_volume ./ (agar + 5);


numPoints = 0;
for i = 1:length(v2)
    numPoints = numPoints + length(v2{i});
end

yMat = zeros(numPoints,1);
A = zeros(numPoints, 8);
ind = 1;

t1 = false;

for i = 1:length(v2)
    kA = nicl_conc(i);
    kB = mncl_conc(i);
    kC = agar(i);
    
    for j = 1:length(v2{i})
       if (t1)
        A(ind,:) = [kA kB kC 0 0 0 1 0 ];
        yMat(ind) = v2{i}(j);
       else
        A(ind,:) = [0 0 0 kA kB kC 0 1];
        yMat(ind) = v2{i}(j);
       end
       ind = ind + 1;
    end
end

figure;
plot(nicl_conc, 1000./R2vals{1},'ro-');
figure
plot(mncl_conc, 1000./R2vals{1},'go-');


A \ yMat

