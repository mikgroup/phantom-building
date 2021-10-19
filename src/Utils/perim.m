function p = perim(boundary)
boundary_dist = 0;
for i = 1:length(boundary)
    if i ~= length(boundary)
        pixel1 = boundary(i,:);
        pixel2 = boundary(i+1,:);
    else
        pixel1 = boundary(i,:);
        pixel2 = boundary(1,:); % when it reaches the last boundary pixel
    end
    pixel_dist = ((pixel1(1,1) - pixel2(1,1)).^2 + (pixel1(1,2) - pixel2(1,2)).^2).^0.5;
    boundary_dist = boundary_dist + pixel_dist;
end
p = boundary_dist;
end