%clear
close all

load T1fit.mat
rng(1);

mask = mask(73:184,:,:,:);
T1est = T1est(73:184,:,:,:);

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
    im = label2rgb(labels(:,:,sl), @jet, [.5 .5 .5]);
    %imshow(label2rgb(labels(:,:,sl), @jet, [.5 .5 .5]))
    st((labels(:,:,sl) > 0).*mean(1000.*T1est(:,:,sl,:),4),[])
    colormap(gca, 'parula')
    colorbar
    axis off
    title('T1 Map (ms)')
    
    bl = 6;
    
    hold on
    for i = 1:length(boundaries{sl})
        boundary = boundaries{sl}{i};
        plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2);
    end
    
    % Find centers of blobs
    stats = regionprops(labels(:,:,sl), 'centroid');
    centers = cat(1,stats.Centroid);
    %plot(centers(:,1), centers(:,2), 'ro');
    
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
        viscircles(centers(blob,:), radiiToUse(blob), 'Color', 'b', 'LineWidth', 1);
        
        circMask = circMask + ...
            double(sqrt((x - centers(blob,1)).^2 + ...
            (y - centers(blob,2)).^2) <= radiiToUse(blob));
        
        row = ceil(centers(blob,1));
        col = ceil(centers(blob,2));
        %h = text(row, col, num2str(labels(col, row, 1)));
        %set(h,'FontSize',14,'FontWeight','bold');
    end
    
    labels(:,:,sl) = labels(:,:,sl) .* circMask;
    hold off 
end
num = max(reshape(labels, [], ns), [], 1).';

%%
slices = [1];
% For 0.11
idx1 = [6, 12, 18, 5, 10, 16, 3, 8, 13, 1, 9, 14, 4, 7, 15, 2, 11, 17];
% For 0.07
%idx1 = [5, 12, 18, 6, 10, 16, 3, 8, 13, 1, 9, 14, 4, 7, 15, 2, 11, 17];

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


%%

niclConc = reshape(repmat(0.25.*[0.7725,5.2487,0.3245,1.3619,8.0507,2.6321]',1,3)',1,18);
mnclConc = reshape(repmat(0.25.*[0.1715,0.0453,0.0986,0.0701,0.3776,0.2463]',1,3)',1,18);
agarWV = repmat([2,2.5,3],1,6)./400;

order = (1:18);
niclConc = niclConc(order);
mnclConc = mnclConc(order);
agarWV = agarWV(order);

numPoints = 0;
for i = 1:length(v2)
    numPoints = numPoints + length(v2{i});
end

yMat = zeros(numPoints,1);
A = zeros(numPoints, 8);
ind = 1;

t1 = true;

for i = 1:length(v2)
    kA = niclConc(i);
    kB = mnclConc(i);
    agarConc = agarWV(i);
    
    for j = 1:length(v2{i})
       if (t1)
        A(ind,:) = [kA kB agarConc 0 0 0 1 0];
        yMat(ind) = v2{i}(j);
       else
        A(ind,:) = [0 0 0 kA kB agarConc 0 1];
        yMat(ind) = v2{i}(j);
       end
       ind = ind + 1;
    end
end



%plot(agarWV, 1000./R1vals{1}, 'ro');

A \ yMat
