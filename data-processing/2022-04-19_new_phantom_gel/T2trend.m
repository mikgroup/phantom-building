close all

load T2fit.mat
rng(10);

mask = T2mask;
mask(T2est > 3) = 0;
mask(T2est < 0.002) = 0;
mask = (flip(mask,1));
T2est = (flip(T2est,1));

labels_cc = zeros(size(mask));
boundaries = cell(length(ns), 1);
SE = strel('diamond',5);
for ii=1:ns
    m0 = mask(:,:,ii);
    m1 = double(imerode(m0, SE));
    m1BW = imbinarize(m1);
    m1BW = bwareaopen(m1BW, 10); % Remove small regions
    
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

%idx1 = [5,10,16,19,2,11,13,20,3,7,17,21,1,8,14,23,4,9,15,22,6,12,18,24];
%idx1 = order_labels(centers);
idx1 = manual_order_labels(centers, labels);
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
arrangement_180 = ([1,7:2:17,21:2:27]);
arrangement_360 = ([4,6:2:18,20:2:27]);
arrangement = arrangement_360;

TPO = [0.4444*ones(1,18),0.8889*ones(1,9)];
acry = [0.4444*ones(1,6),0.8889*ones(1,6),1.3333*ones(1,6),0.4444*ones(1,6),0.8889*ones(1,3)];
PEGDA = [repmat([0.4444,0.4444,0.8889,0.8889,1.3333,1.3333],[1,4]),[0.4444,0.4444,0.8889]];
water = [2.6667,2.6667, 2.2222, 2.2222, 1.7778,1.7778, 2.2222,2.2222,1.7778,1.7778,1.3333,...
    1.3333,1.7778,1.7778,1.3333,1.3333,0.8889,0.8889,2.2222,2.2222,1.7778,1.7778,1.3333,1.3333,1.7778,1.7778,1.3333];

TPO = TPO(arrangement);
acry = acry(arrangement);
PEGDA = PEGDA(arrangement);
water = water(arrangement);

TPO_conc = TPO./(TPO+acry+PEGDA+water);
acry_conc = acry./(TPO+acry+PEGDA+water);
PEGDA_conc = PEGDA./(TPO+acry+PEGDA+water);


numPoints = 0;
for i = 1:length(arrangement)
    numPoints = numPoints + length(v2{arrangement(i)});
end

yMat = zeros(numPoints,1);
A = zeros(numPoints, 8);
ind = 1;

t1 = false;

for i = length(arrangement)
    kA = TPO_conc(i);
    kB = acry_conc(i);
    kC = PEGDA_conc(i);
    
    for j = 1:length(v2{arrangement(i)})
       if (t1)
        A(ind,:) = [kA kB kC 0 0 0 1 0 ];
        yMat(ind) = v2{arrangement(i)}(j);
       else
        A(ind,:) = [0 0 0 kA kB kC 0 1];
        yMat(ind) = v2{arrangement(i)}(j);
       end
       ind = ind + 1;
    end
end


figure
plot(TPO_conc, R2vals{1}(arrangement),'r*');
hold on
plot(PEGDA_conc, R2vals{1}(arrangement),'b*');
plot(acry_conc, R2vals{1}(arrangement),'k*')
hold off
%A \ yMat
%%


function idx = manual_order_labels(centers, labels)

figure, imshow(labels.^(1/2),[])

idx = zeros(size(centers,1),1);

n=1;
while(1)
    [x,y] = ginput(1);
    dist = sqrt(sum((centers-[x,y]).^2,2));
    tmp = find(dist==min(dist));
    if sum(tmp==idx)==0
        text(centers(tmp,1),centers(tmp,2),sprintf('%d',n),'Color','red', 'FontSize',12);
        idx(n) = tmp;
        n=n+1;
        if n>size(centers,1)
            break;
        end
    end
end

end
    
    
