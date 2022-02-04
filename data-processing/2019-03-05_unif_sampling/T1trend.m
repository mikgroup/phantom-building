%clear
close all

load T1fit.mat
rng(1);

mask(T1est > 2.8) = 0;
mask(T1est < 0.100) = 0;
mask = flip(flip(mask,1),2);

T1est = flip(flip(T1est,1),2);

labels_cc = zeros(size(mask));
boundaries = cell(length(ns), 1);
SE = strel('diamond',2);
for ii=1:ns
    m0 = mask(:,:,ii);
    m1 = double(imerode(m0, SE));
    m1BW = imbinarize(m1);
    m1BW = bwareaopen(m1BW, 20); % Remove small regions
    
    % Labels and creates cell array of boundaries
    [B, L] = bwboundaries(m1BW); 
    labels_cc(:,:,ii) = L;
    boundaries{ii} = B;
end

labels = labels_cc;
clear labels_cc m0 m1 L B SE

%%

for sl = 1:size(labels,3)
    % Plots boundaries
    % make sure these are not touching
    figure
    
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

idx1 = order_labels(centers);
idxs = {idx1};

R1vals = cell(1, length(slices));

map = T1est;
stds = zeros(1,length(idxs));
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
        stds(jj) = std(1000./x3);
        %hist(x3,50)
    end
    
    R1vals{ii} = v;
end

%%

%%

nicl_sol_volume = [1.8279,0.2359,2.7895,1.1975,0.4127,3.2427,1.6507,0.8659,0.3985,0.0883,1.9144,1.1295,0.6621,0.352,0.1312,1.302,0.8346,0.5245,0.3036,0.1384,0.9562,0.6461,0.4253,0.26,0.1317,0.3504,0.2221,0.1196];
mncl_sol_volume = [0.5625,0.6246,0.2381,0.3002,0.3308,0.0852,0.1473,0.1779,0.1961,0.2082,0.0583,0.0889,0.1071,0.1192,0.1279,0.0307,0.049,0.0611,0.0697,0.0761,0.0079,0.02,0.0286,0.0351,0.0401,0.0046,0.0096,0.0136];
agar = [16.0689, 16.0675, 16.1007, 16.0766, 16.1011, 16.1116, 16.0564, 16.096, 16.0982, 16.1234, 16.0925, 16.0991, 16.1088, 16.0938,...
    16.0765, 16.1229, 16.0906, 16.0833, 16.1425, 16.0867, 16.1055, 16.2060, 16.0745, 16.1096, 16.0821, 16.1435, 16.1514, 16.1533];

vial_placement = [3, 27, 14, 6, 26, 7, 22, 10, 23, 18, 8, 17, 9, 1, 19, 20, 15, 13, 4, 25, 16, 5, 11, 2, 21, 12, 24, 28];
nicl_sol_volume = nicl_sol_volume(vial_placement);
mncl_sol_volume = mncl_sol_volume(vial_placement);
agar = agar(vial_placement);


nicl_conc = nicl_sol_volume ./ (agar + 5);

mncl_conc = mncl_sol_volume ./ (agar + 5);

numPoints = 0;
for i = 1:length(v2)
    numPoints = numPoints + length(v2{i});
end

yMat = zeros(numPoints,1);
A = zeros(numPoints, 8);
ind = 1;

t1 = true;

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



plot(nicl_conc, R1vals{1},'ro-');
hold on
plot(mncl_conc, R1vals{1},'ro-');
hold off


%plot(agarWV, 1000./R1vals{1}, 'ro');

%A \ yMat
