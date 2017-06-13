load T1fit.mat
rng(1);

T1est_orig = T1est;
T1est = flip(flip(T1est, 1), 2);
mask_orig = mask;
mask = flip(flip(mask, 1), 2);

mask = mask;

labels_cc = zeros(size(mask));
boundaries = cell(length(ns), 1);
SE = strel('diamond',2);
for ii=1:ns
    m0 = mask(:,:,ii);
    m1 = double(imerode(m0, SE));
    m1BW = imbinarize(m1);
    
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
    imshow(label2rgb(labels(:,:,sl), @jet, [.5 .5 .5]))
    
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
idx1 = [6,5,3,4,1,2,10,7,8,9,11,12];

idxs = {idx1};

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
        i2 = find(x2 <= 1.99*median(x2), 1, 'last');
        x3 = x2(i1:i2); % throw out outliers
        v(jj) = mean(x3);
        v2{jj} = x3;
    end
    
    R1vals{ii} = v;
end
hold off


%%

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

figure;
plot(vialWV, R1vals{1});

%A \ yMat
