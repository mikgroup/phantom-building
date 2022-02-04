% JT 07-2016
%clear;
%close all;
T1imgs_raw = sqreadcfl('b1_cimg');
T1imgs_raw = permute(T1imgs_raw, [1, 2, 5, 3, 4]);
[ny, nz, ns, nc, nt] = size(T1imgs_raw);
%%

B1imgs = squeeze(mean(T1imgs_raw,4));
B1est = acos(B1imgs(:,:,2)./(2*B1imgs(:,:,1)));