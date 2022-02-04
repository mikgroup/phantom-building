%clear
close all

load T1fit.mat
rng(1);

%T1est = mean(T1est, 4);
mask(T1est > 3.2) = 0;
mask(T1est < 0.01) = 0;

% imshow(T1est(:,:,1,1))
% r = drawfreehand;
% mask = mask & createMask(r);

% imshow(mask(:,:,1,1))
% save('T1_mask.mat','mask');

load('T2_mask.mat')

mask = (flip(mask,1));

T1est = (flip(T1est,1));

labels_cc = zeros(size(mask));
boundaries = cell(length(ns), 1);
SE = strel('diamond',3);
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

%idx1 = [6    12     1     7    14     2     9    15     3    10    16    11     5    13     8];
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

arrangement = ["C4", "M2", "C2", "C1", "M1", "M4", "M2", "M4", "M3", "L1.5", "L3", "M3", "C3","L1","M1","C3","C2","M1","C4","L2.5","C0","C0","L3.5","L4","C2","L2","L0.5","C1","C0","M3","C4","M2","C3","C1","M4"];

ion = "M";

mncl_pos = arrangement.contains(ion) | arrangement.contains("C0");
mncl_c = erase(erase(arrangement(mncl_pos),"M"),"C");
mncl_conc = zeros(size(mncl_c));
for i = 1:size(mncl_c,2)
    mncl_conc(i) = str2num(mncl_c(i));
end

mncl_r1 = R1vals{1}(mncl_pos);


ion = "C";

cono_pos = arrangement.contains(ion);
cono_c = erase(arrangement(cono_pos),ion);
cono_conc = zeros(size(cono_c));
for i = 1:size(mncl_c,2)
    cono_conc(i) = str2num(cono_c(i));
end

cono_r1 = R1vals{1}(cono_pos);


figure;
plot(mncl_conc, mncl_r1,'o','MarkerSize', 7,'LineWidth',2);

p = polyfit(mncl_conc', mncl_r1,1);
x = linspace(0,max(mncl_conc));
y = polyval(p, x);

hold on
plot(x,y,'LineWidth',2)
plot(cono_conc, cono_r1,'o','MarkerSize', 7,'LineWidth',2);
hold off

xlabel('NiCl_2 Gel Volume (mL)')
ylabel('R1 (1/s)')
faxis(gca, 20);

figure;
plot(mncl_conc, 1000./mncl_r1,'ro','MarkerSize', 7,'LineWidth',2);
hold on
plot(x,1000./y,'LineWidth',2)
hold off

xlabel('NiCl_2 Gel Volume (mL)')
ylabel('T1 (ms)')
faxis(gca, 20);

%% Cure time

cure_pos = arrangement.contains("L");
cure_t = erase(arrangement(cure_pos),"L");
cure_time = zeros(size(cure_t));
for i = 1:size(cure_t,2)
    cure_time(i) = str2num(cure_t(i));
end

cure_r1 = R1vals{1}(cure_pos);

figure;
plot(cure_time, 1000./cure_r1,'o','MarkerSize', 7,'LineWidth',2);
xlabel('Curing Time (min)')
ylabel('T1 (ms)')
title('T1 Recovery vs Curing Time')
faxis(gca, 20);

%A \ yMat
