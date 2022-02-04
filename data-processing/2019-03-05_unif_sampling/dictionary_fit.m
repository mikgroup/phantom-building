function [T2, M0, T1, B1, idx] = dictionary_fit(y, D, dims, T2_arr, T1_arr, B1_arr)
zz = D*y;
[~, idx] = max(abs(zz));
M0 = zz(idx);
[i1, i2, i3] = ind2sub(dims, idx);
T2 = T2_arr(i3);
T1 = T1_arr(i2);
B1 = B1_arr(i1);

end

