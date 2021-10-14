clear
close all
load 'brain.mat'

% Slices for each dimension
%1: 84 
%2: 103
%3: 94 (84) (64)

slice = 64;
scale = 5;

brainSlice = flip(imrotate(squeeze(brain(:,:,94)),90)')';
st(round(imresize(brainSlice,scale)), []); ...
    %colormap(gca, 'parula')       
brainSlice(brainSlice == 8) = 2; % Assign glial cells to CSF
brainSlice = brainSlice .* (brainSlice <= 3);
brainSliceScale = imresize(brainSlice, scale);
brainSliceScale = round(brainSliceScale);
%%
mask = brainSlice > 0;
mask = round(imresize(mask, scale));
brainSliceScale = brainSliceScale .* mask;


st(brainSliceScale.* (brainSliceScale >= 1 & brainSliceScale <= 3), []);...
    colormap(gca, 'parula')

brainSliceScale = brainSliceScale'';
st(brainSliceScale.* (brainSliceScale >= 1 & brainSliceScale <= 3), []);...
    colormap(gca, 'parula')

brainSliceScale = brainSliceScale'';

%figure;
si = size(brainSliceScale);
axis([0 si(1) 0 si(2)])
axis equal
axis off
set(gcf,'color','w');
%imshow(zeros(size(brainSliceScale)));
hold on
%colors = ['r', 'b', 'g'];
colors = ['k', 'k', 'k'];
bounds = ones(size(brainSliceScale)).*0;
ind = 1;
perims = [];
allBounds = {};
allBoundsInd = 1;
cents = [];

hold on
for reg = 1:3  
    [B, L] = bwboundaries(brainSliceScale == reg, 'noholes');
    c = regionprops(brainSliceScale == reg, 'centroid');
    cents = [cents; cat(1, c.Centroid)];
    for i = 1:length(B)
        %B{i} = [B{i}; B{i}(2,:)]; % add point to avoid discontinuity in plot
        p = perim(B{i});
        perims(ind) = p;
        ind = ind + 1;
        if p > (100 * scale) % change this number to remove small perimeters
            allBounds{allBoundsInd} = B{i};
            allBoundsInd = allBoundsInd + 1;
            plot(B{i}(:,2), B{i}(:,1),colors(reg),'LineWidth', 1.5);
            for j = 1:length(B{i})
                bounds(B{i}(j,1),B{i}(j,2)) = 255;%bounds(B{i}(j,1),B{i}(j,2)) + 1;
            end
        end
    end
end
hold off
st(bounds,[]), colormap(gca, 'parula')
xsum = sum(bounds);
x_min = find(xsum > 0, 1, 'first');
x_max = find(xsum > 0, 1, 'last');
ysum = sum(bounds, 2);
y_min = find(ysum > 0, 1, 'first');
y_max = find(ysum > 0, 1, 'last');
%%
conv_to_mm = 1/scale;%183/(x_max - x_min);
allBoundsScaled = cell(size(allBounds));
for i = 1:length(allBounds)
    allBoundsScaled{i} = allBounds{i}.*conv_to_mm;
end

b = allBoundsScaled{1};
[x,y] = offsetCurve(b(:,1), b(:,2), 1);
d = [x',y'];
d(end+1,:) = d(1,:);
allBoundsScaled{1} = d;

com_x = sum(d(:,2))/length(d);
com_y = sum(d(:,1))/length(d);

targ_x = 100;%210 - 95.1671;
targ_y = 100;%157.5 - 56.0216;

% for i = 1:length(allBoundsScaled)
%     allBoundsScaled{i} = flip(allBoundsScaled{i}, 2); %flip x,y
%     allBoundsScaled{i}(:,1) = allBoundsScaled{i}(:,1) + targ_x - com_x;
%     allBoundsScaled{i}(:,2) = allBoundsScaled{i}(:,2) + targ_y - com_y;
% end

figure
hold on
%axis square
%axis([0 420 0 315])
for i = 1:length(allBoundsScaled)
    allBoundsScaled{i}(end + 1,:) = allBoundsScaled{i}(1,:);
    plot(allBoundsScaled{i}(:,2)', -allBoundsScaled{i}(:,1)', 'k', 'linewidth', 1.5);
end
axis equal off
%gen_gcode(allBoundsScaled,'phantom_aspect.gcode')


