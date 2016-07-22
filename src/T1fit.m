% JT 07-2016

imgs_raw = sqreadcfl('~/Desktop/2016-07-19_phantom-test/sandbox/t1_tir_img');
[ny, nz, ns, T, R] = size(imgs_raw);
imgs_raw = squeeze(imgs_raw(:,:,:,1,:)); % only need first spin echo
%%
d1 = squeeze(dimnorm(imgs_raw(:,:,3,:), 4));
mask = d1 > .05*max(d1(:));

inv_times = [250, 499, 1000] * 1e-3;

proton1 = zeros(ny, nz, ns);
proton2 = zeros(ny, nz, ns);
T1est = zeros(ny, nz, ns);

%% 15 sec on [320, 322, 5, 3] using Joelle's awesome method
% S(TIn) = (a + b exp(-TIn/T1)), where a and b are complex-valued
myfun0 = @(v, inv_times) (v(1)+1j*v(2)) + (v(3)+1j*v(4))*exp(-inv_times*v(5));
myfun = @(v, inv_times) [real(myfun0(v, inv_times)), imag(myfun0(v, inv_times))];

extra.T1Vec = (1:5000)*1e-3;
extra.tVec = inv_times;
nlsS = getNLSStruct(extra);

tic
s = 60;
p = s/(ny*nz*ns);
fprintf(1,'|%s|\n|\n',repmat('-',1,s));
parfor ii=1:ny*nz*ns
    if rand < p
        fprintf(1,'\b.\n'); % \b is backspace
    end
    [yy, zz, cc, ss] = ind2sub([ny, nz, ns], ii);
    mm = mask(yy, zz);
    if mm == 0
        proton1(ii) = 0;
        proton2(ii) = 0;
        T1est(ii) = 0;
    else
        y_cplx = squeeze(imgs_raw(yy,zz,ss,:));


        [T1, v1, v2, res] = rdNls(y_cplx, nlsS);
        
        proton1(ii) = v1;
        proton2(ii) = v2;
        T1est(ii) = T1;
    end
end
toc
%%
figure(1); hist(T1est(T1est~=0)*1000, 100); faxis
st(1000*bsxfun(@times, mean(mean(T1est,3),4), mask), [0, 1000]); colormap('default'), colorbar;
title('T1 map (ms)'); faxis
st(1000*bsxfun(@times, T1est, mask), [0, 1000]); colormap('default'), colorbar;
title('T1 map (ms)'); faxis
stc(bsxfun(@times, proton1, mask)); colormap('default'), colorbar;
title('mag1 map'); faxis
stc(bsxfun(@times, proton2, mask)); colormap('default'), colorbar;
title('mag2 map'); faxis
