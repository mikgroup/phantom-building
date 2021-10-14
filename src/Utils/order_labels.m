function res = order_labels(centers)

    % Sort the y center positions
    [centers_y_sort, inds] = sort(centers(:,2));
    edges = diff(centers_y_sort, 1) > mean(diff(centers_y_sort), 1);
    
    
    rows = [0; find(edges); size(centers,1)];
    
    curr_ind = 1;
    centers_x_sort = centers(:,1);
    centers_x_sort = centers_x_sort(inds);
    
    sub_sort_inds = [];
    for i = 2:size(rows)
        [~, x_inds] = sort(centers_x_sort(rows(i - 1)+1:rows(i)));
        
        sub_sort_inds = [sub_sort_inds; x_inds+ rows(i - 1)];
        curr_ind = rows(i);
    end
    res = inds(sub_sort_inds)
end