%clear
close all

load T1fit.mat
rng(1);

mask = flip(flip(mask,1),2);
T1est = flip(flip(T1est,1),2);

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
% For 0.11
idx1 = [4,5,3,2,1,8,7,6,9];
% For 0.07
%idx1 = [5, 12, 18, 6, 10, 16, 3, 8, 13, 1, 9, 14, 4, 7, 15, 2, 11, 17];

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

nicl_conc = [0.5,1,1.5,0.5,1,1.5,0.5,1,1.5] .* (5/20);
mncl_conc = [0.25,0.25,0.25,0.5,0.5,0.5,0.75,0.75,0.75] .* (5/20);


%R1vals{1} = sort(R1vals{1});
numPoints = 0;
for i = 1:length(v2)
    numPoints = numPoints + length(v2{i});
end

figure;
hold on
t1 = 1000./R1vals{1};
for i = 1:3
    inds = (1+3*(i-1)):(3+3*(i-1));
    t1_cur = t1(inds);
    plot(nicl_conc(inds), 1./t1_cur, 'o-');
    a = polyfit(nicl_conc(inds), 1./t1_cur',1);
    x_vals = min(nicl_conc):0.01:max(nicl_conc);
    y_vals = a(1).*x_vals + a(2);
    plot(x_vals, y_vals);
    a(1)
end
hold off

yMat = zeros(numPoints,1);
A = zeros(numPoints, 6);
ind = 1;

t1 = true;

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


%plot(nicl_conc, R1vals{1},'ro-');


%plot(agarWV, 1000./R1vals{1}, 'ro');

%A \ yMat
