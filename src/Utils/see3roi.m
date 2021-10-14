function [ a ] = see3roi(IM, noabs, rescale, figpointer, count)

use_abs = false;
if nargin < 2 || noabs == 0
    use_abs = true;
end

if nargin < 3
    rescale = false;
end

if nargin < 4
    figpointer = 22;
end

if nargin < 5
    count = inf
end

button=1;
im = sqrt(sum(abs(IM.^2),3));


if rescale
    imsort = sort(abs(IM(:)), 'ascend');
    ymin = imsort(1);
    ymax = imsort(round(end - end/1000));
end

figmain = figpointer;
figpointer = figpointer + 1;

while count > 0
    count = count - 1;
    figure(figmain), imshow(im,[]), impixelinfo

    h = imellipse;
    BW = createMask(h);
    
    a = 0 * IM(1,1,:);
    
    figure(11);
    for ii=1:size(IM,3)
        IM1 = IM(:,:,ii);
        if use_abs
            a(ii) = mean(abs(IM1(BW)));
        else
            a(ii) = mean(IM1(BW));
        end
    end
    
    a = a(:);
    
    figure(figpointer), plot(a);
    if rescale
        ylim([ymin, ymax]);
    end
%     title(sprintf('(%d, %d)',round(y),round(x)));
    if button == 3
        figpointer = figpointer + 1;
    end
    close(11);
end

