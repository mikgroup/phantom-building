function [T2, M0, idx] = dictionary_fit_T2_only(y, echo_train_modulation, T2_arr, T1_arr, B1_arr, T1, B1)
[~, idx_T1] = min(abs(T1 - T1_arr));
[~, idx_B1] = min(abs(B1 - B1_arr));

ETL = size(echo_train_modulation, 5);
D = reshape(echo_train_modulation(:,idx_B1,idx_T1,:,:), [], ETL); D = bsxfun(@rdivide, D, dimnorm(D,2));
zz = D*y;
[~, idx] = max(abs(zz));
M0 = zz(idx);
T2 = T2_arr(idx);
end
