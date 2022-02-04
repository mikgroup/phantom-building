% JT 07-2016
clear
load dictionary_data.mat
T2imgs_raw = sqreadcfl('se_mc_cimg');
T2imgs_raw = permute(T2imgs_raw, [1, 2, 5, 3, 4]); % only map specific slice
[ny, nz, ns, nc, nt] = size(T2imgs_raw);
%%
d1 = dimnorm(T2imgs_raw, 5);
%d2 = dimnorm(d1, 4);
d2 = d1;
T2mask = d2 > .01*max(d2(:));

echo_times = [8.8:8.8:32*8.8]*1e-3; %
echo_times = echo_times(2:2:2*16);

%%
proton = zeros(ny, nz, ns, nc);
T2est = zeros(ny, nz, ns, nc);

%% 300 sec on [256, 256, 6, 32]
myfun0 = @(v, echo_times) complex(v(1), v(2)) * exp(-echo_times*v(3));
myfun = @(v, echo_times) [real(myfun0(v, echo_times)), imag(myfun0(v, echo_times))];
ops = optimoptions('lsqcurvefit','TolFun', 1e-10, 'Display', 'off');
LB = [-inf, -inf, 1/10];
UB = [inf, inf, 1/5e-3];
tic
s = 60;
p = s/(ny*nz*ns*nc);
rvals = rand(ny*nz*ns*nc, 1);
fprintf(1,'|%s|\n|\n',repmat('-',1,s));
parfor ii=1:ny*nz*ns*nc
    if rvals(ii) < p
        fprintf(1,'\b.\n'); % \b is backspace
    end
    [yy, zz, ss, cc] = ind2sub([ny, nz, ns, nc], ii);
    mm = T2mask(yy, zz, ss);
    if mm == 0
        proton(ii) = 0;
        T2est(ii) = 0;
    else
        % -- complex-valued fit  -- %
        y_cplx = squeeze(T2imgs_raw(yy,zz,ss,cc,:)); 
        
        [T2_curr, M0, ~, ~, ~] = dictionary_fit(y_cplx, D, dims, T2_tse_arr, T1_tse_arr, B1_scaling_arr);  
        
        proton(ii) = M0;
        T2est(ii) = T2_curr;
        
        % -- magnitude fit -- %
        %             y = abs(y_cplx);
        %             F = fit(echo_times.', y, 'exp1');
        %             x = [F.b, F.a];
        %             proton(yy, zz) = F.a;
        %             T2est(yy, zz) = -1/F.b;
    end
    %     toc
end
toc

%%
T2vals = T2est(T2est~=0);
figure(1); hist(T2vals(T2vals<.2)*1000, 100); faxis
%
st(1000*bsxfun(@times, T2est, T2mask), [0, 160]); colormap('parula'), colorbar;
title('T2 map (ms)'); faxis
stc(bsxfun(@times, proton, T2mask)); colormap('parula'), colorbar;
title('proton density'); faxis

save T2fit_mc_noam.mat
