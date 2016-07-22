% JT 07-2016
imgs_raw = sqreadcfl('se_mc_img');
[ny, nz, ns, T] = size(imgs_raw);
%%

mask = dimnorm(squeeze(imgs_raw(:,:,3,:)),3) > .05*max(reshape(dimnorm(imgs_raw(:,:,3,:),3),1,[]));
% mask = ones(size(imgs(:,:,1)));

TE = 13.6e-3;
echo_times_full = TE:TE:(TE+(T-1)*TE);
echo_times = echo_times_full(2:end); % throw out first echo

proton = zeros(ny, nz, ns);
T2est = zeros(ny, nz, ns);

%% 300 sec on [256, 256, 6, 32]
myfun0 = @(v, echo_times) complex(v(1), v(2)) * exp(-echo_times*v(3));
myfun = @(v, echo_times) [real(myfun0(v, echo_times)), imag(myfun0(v, echo_times))];
ops = optimoptions('lsqcurvefit','FunctionTolerance', 1e-10, 'Display', 'off');
LB = [-inf, -inf, 0];
UB = [];
tic
s = 60;
p = s/(ny*nz*ns);
fprintf(1,'|%s|\n|\n',repmat('-',1,s));
parfor ii=1:ny*nz*ns
    if rand < p
        fprintf(1,'\b.\n'); % \b is backspace
    end
    [yy, zz, ss] = ind2sub([ny, nz, ns], ii);
    mm = mask(yy, zz);
    if mm == 0
        proton(ii) = 0;
        T2est(ii) = 0;
    else
        % -- complex-valued fit  -- %
        y_cplx = squeeze(imgs_raw(yy,zz,ss,2:end));
        y = [real(y_cplx), imag(y_cplx)];
        v0 = [1e-4, 1e-4, 1/60e-3];
        [v, resnorm, res, exitflag, output] = ...
            lsqcurvefit(myfun, v0, echo_times.', y, LB, UB, ops);
        
        proton(ii) = complex(v(1), v(2));
        T2est(ii) = 1/v(3);
        
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
% figure(1); hist(T2est(T2est~=0)*1000, 100); faxis
st(1000*bsxfun(@times, T2est, mask), [0, 100]); colormap('default'), colorbar;
title('T2 map (ms)'); faxis
stc(bsxfun(@times, proton, mask)); colormap('default'), colorbar;
title('proton density'); faxis
