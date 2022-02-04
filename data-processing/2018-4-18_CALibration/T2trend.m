%clear
close all

load T2fit.mat
rng(10);

mask = flip(flip(mask,1),2);
T2est = flip(flip(T2est,1),2);

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

idx1 = [6,5,4,2,1,3,10,12,11,9,7,8,14,16,15,13];
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

nicl_conc = [0.286581068322652,1.08682719060280,2.88240062254969,5.46372628537846,3.96710074968355,0.338944568073161,1.02280670187358,2.58745579946568,0.980055598326268,0.840704376575678,1.48936864670864,0.225736299854347,0.604267823935120,0.997940843108945,0.726363973687932,0.566535024002472];
mncl_conc = [0.102889930841530,0.126251274743614,0.103113390417274,0.0434947705580906,0.0274454171902448,0.0534401141235837,0.0439943334967211,0.0250527706659723,0.0325778686107493,0.0264027331107057,0.00861998392970116,0.0299801773122196,0.0184406387649377,0.00830499125908914,0.00576205584772758,0.00398619275306566];

numPoints = 0;
for i = 1:length(v2)
    numPoints = numPoints + length(v2{i});
end

yMat = zeros(numPoints,1);
A = zeros(numPoints, 6);
ind = 1;

t1 = false;

for i = 1:length(v2)
    kA = nicl_conc(i);
    kB = mncl_conc(i);
    
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

figure;
plot(nicl_conc, R2vals{1},'ro-');
hold on
plot(mncl_conc, R2vals{1},'go-');
hold off

A \ yMat
