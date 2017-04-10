% JT 07-2016

T1imgs_raw = sqreadcfl('ir_cimg');
T1imgs_raw = permute(T1imgs_raw, [1, 2, 5, 3, 4]);
[ny, nz, ns, nc, nt] = size(T1imgs_raw);
%%
d1 = dimnorm(T1imgs_raw, 5);
% d2 = dimnorm(d1, 4);
d2 = d1;
mask = d2 > .11*max(d2(:));
st(mask)
% clear d1 d2
% inv_times = [250, 499, 1000] * 1e-3; % 2016-07-19_phantom-test
inv_times = [100, 500, 900, 1300, 2000] * 1e-3; % 2016-08-01_phantom-test

%%
proton1 = zeros(ny, nz, ns, nc);
proton2 = zeros(ny, nz, ns, nc);
T1est = zeros(ny, nz, ns, nc);

%% 15 sec on [320, 322, 5, 3] using Joelle's awesome method
% S(TIn) = (a + b exp(-TIn/T1)), where a and b are complex-valued
myfun0 = @(v, inv_times) (v(1)+1j*v(2)) + (v(3)+1j*v(4))*exp(-inv_times*v(5));
myfun = @(v, inv_times) [real(myfun0(v, inv_times)), imag(myfun0(v, inv_times))];

extra.T1Vec = (1:5000)*1e-3;
extra.tVec = inv_times;
nlsS = getNLSStruct(extra);

tic
s = 60;
p = s/(ny*nz*ns*nc);
fprintf(1,'|%s|\n|\n',repmat('-',1,s));
parfor ii=1:ny*nz*ns*nc
    if rand < p
        fprintf(1,'\b.\n'); % \b is backspace
    end
    [yy, zz, ss, cc] = ind2sub([ny, nz, ns, nc], ii);
    mm = mask(yy, zz, ss);
    if mm == 0
        proton1(ii) = 0;
        proton2(ii) = 0;
        T1est(ii) = 0;
    else
        y_cplx = squeeze(T1imgs_raw(yy,zz,ss,cc,:));


        [T1, v1, v2, res] = rdNls(y_cplx, nlsS);
        
        proton1(ii) = v1;
        proton2(ii) = v2;
        T1est(ii) = T1;
    end
end
toc
%%
figure(1); hist(T1est(T1est~=0)*1000, 100); faxis
st(1000*bsxfun(@times, mean(mean(T1est,3),4), mask), [0, 1000]); colormap('parula'), colorbar;
title('T1 map (ms)'); faxis
%%
st(1000*bsxfun(@times, T1est, mask), [0, 1000]); colormap('parula'), colorbar;
title('T1 map (ms)'); faxis
stc(bsxfun(@times, proton1, mask)); colormap('parula'), colorbar;
title('mag1 map'); faxis
stc(bsxfun(@times, proton2, mask)); colormap('parula'), colorbar;
title('mag2 map'); faxis
