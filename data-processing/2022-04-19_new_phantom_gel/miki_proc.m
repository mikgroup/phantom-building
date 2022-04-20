%%

addpath ~/'Google Drive'/'Shared drives'/mikGroupDrive/Home/mlustig/Matlab/ESPIRiT/
addpath ~/'Google Drive'/'Shared drives'/mikGroupDrive/Home/mlustig/Matlab/ESPIRiT/utils/
addpath ~/'Google Drive'/'Shared drives'/mikGroupDrive/Home/mlustig/Matlab/ESPIRiT/ESPIRiT_code/
addpath ~/'Google Drive'/'Shared drives'/mikGroupDrive/Home/mlustig/Matlab/orchestra-sdk-2.0-1-2.matlab
addpath ~/'Google Drive'/'Shared drives'/mikGroupDrive/Home/mlustig/Matlab/orchestra-sdk-2.0-1-2.matlab/Scripts/
addpath ~/'Google Drive'/'Shared drives'/mikGroupDrive/Home/mlustig/Matlab/orchestra-sdk-2.0-1-2.matlab/Examples/
addpath ~/'Google Drive'/'Shared drives'/mikGroupDrive/Home/mlustig/Matlab/phantom-building/src/Utils

%%
% file names
files = dir('P*7');

%%
%load k-space and extract TI, TE -- match name to header. 

c_TE = 0;
c_TI = 0;
kspaces_TI=[];
kspaces_TE=[];
for n=1:length(files)
    pfile = GERecon('Pfile.Load',files(n).name,'No-Anonymize');
    header = GERecon('Pfile.Header',pfile);
    if strfind(header.SeriesData.se_desc,'TI')
        c_TI = c_TI+1;
        nameTI = str2num(header.SeriesData.se_desc(strfind(header.SeriesData.se_desc,"=")+1:end));
        TI(c_TI) = floor(header.ImageData.ti/1000);
        if nameTI~= TI(c_TI);
            disp(sprintf('Warning! name TI=%d does not match sequence TI=%d',nameTI,TI(c_TI)));
        else
            disp(sprintf('%s:TI=%d',files(n).name,TI(c_TI)));
        end
        kspaces_TI(:,:,:,c_TI) = pfile2kspace(pfile);
    end
    if strfind(header.SeriesData.se_desc,'TE')
        c_TE = c_TE+1;
        TE(c_TE) = floor(header.ImageData.te/1000);
        nameTE = str2num(header.SeriesData.se_desc(strfind(header.SeriesData.se_desc,"=")+1:end));
        if nameTE~= TE(c_TE);
            disp(sprintf('Warning name TE=%d does not match sequence TE=%d',nameTE,TE(c_TE)));
        else
            disp(sprintf('%s: TE=%d',files(n).name,TE(c_TE)));
        end
        kspaces_TE(:,:,:,c_TE) = pfile2kspace(pfile);
    end
end

    
%%
% Sort TE and TI

[TI,idx] = sort(TI);
kspaces_TI = kspaces_TI(:,:,:,idx);

[TE,idx] = sort(TE);
kspaces_TE = kspaces_TE(:,:,:,idx);

%%
% compute sensitivity maps 

% choose scan to compute maps for all data (usually, first TI)
DATA = kspaces_TI(:,:,:,1);
im = ifft2c(DATA);

[sx,sy,Nc] = size(DATA);
ncalib = 24; % use 24 calibration lines to compute compression
ksize = [6,6]; % kernel size
% Threshold for picking singular vercors of the calibration matrix
% (relative to largest singlular value.

eigThresh_1 = 0.02;

% threshold of eigen vector decomposition in image space. 
eigThresh_2 = 0.95;

% crop a calibration area
calib = crop(DATA,[ncalib,ncalib,Nc]);

% compute Calibration matrix, perform 1st SVD and convert singular vectors
% into k-space kernels

[k,S] = dat2Kernel(calib,ksize);
idx = max(find(S >= S(1)*eigThresh_1));
 
% crop kernels and compute eigen-value decomposition in image space to get
% maps
[M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);

ESP = ESPIRiT(M(:,:,:,end));
P = ESP'*im;
err = sos(im-ESP*P);

figure, imshow(cat(2,abs(P),abs(err)).^(1/3),[])

%%
% coil combine image

images_TI = [];
for n=1:length(TI)
    images_TI(:,:,1,1,n) = ESP'*(ifft2c(kspaces_TI(:,:,:,n)));
end

images_TE = [];
for n=1:length(TE)
    images_TE(:,:,1,1,n) = ESP'*(ifft2c(kspaces_TE(:,:,:,n)));
end
    

%% 
% Gradient warp params

corners = GERecon('Pfile.Corners', 1,pfile);
orientation = GERecon('Pfile.Orientation', 1,pfile);
images_TI_gw = [];
for n=1:length(TI)
    tmp = ifft2c(zpad(fft2c(images_TI(:,:,1,1,n)),[512,512]));
    tmp = GERecon('Gradwarp', real(tmp), corners,'XRMW')+1j*GERecon('Gradwarp', imag(tmp), corners,'XRMW');
    tmp = crop(ifft2c(crop(fft2c(tmp),[size(images_TI,1),size(images_TI,1)])),[size(images_TI,1),size(images_TI,2)]);
    images_TI_gw(:,:,1,1,n) = tmp;
end



images_TE_gw = [];
for n=1:length(TE)
    tmp = ifft2c(zpad(fft2c(images_TE(:,:,1,1,n)),[512,512]));
    tmp = GERecon('Gradwarp', real(tmp), corners,'XRMW')+1j*GERecon('Gradwarp', imag(tmp), corners,'XRMW');
    tmp = crop(ifft2c(crop(fft2c(tmp),[size(images_TE,1),size(images_TE,1)])),[size(images_TE,1),size(images_TE,2)]);
    images_TE_gw(:,:,1,1,n) = tmp;
end
%%
images_TI_gw = images_TI;
images_TE_gw = images_TE;


%% 
% T1 fit:

[ny, nz, ns, nc, nt] = size(images_TI_gw);
d1 = dimnorm(images_TI_gw, 5);
d2 = d1;
mask = d2 > .01*max(d2(:));
inv_times = TI* 1e-3;

proton1 = zeros(ny, nz, ns, nc);
proton2 = zeros(ny, nz, ns, nc);
T1est = zeros(ny, nz, ns, nc);

myfun0 = @(v, inv_times) (v(1)+1j*v(2)) + (v(3)+1j*v(4))*exp(-inv_times*v(5));
myfun = @(v, inv_times) [real(myfun0(v, inv_times)), imag(myfun0(v, inv_times))];

extra.T1Vec = (1:10000)*1e-3;
extra.tVec = inv_times;
nlsS = getNLSStruct(extra);

tic
s = 60;
p = s/(ny*nz*ns*nc);
fprintf(1,'|%s|\n|\n',repmat('-',1,s));
for ii=1:ny*nz*ns*nc
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
        y_cplx = squeeze(images_TI_gw(yy,zz,ss,cc,:));


        [T1, v1, v2, res] = rdNls(y_cplx, nlsS);
        
        proton1(ii) = v1;
        proton2(ii) = v2;
        T1est(ii) = T1;
    end
end
toc

save T1fit.mat
%%
% Display T1

figure(1); hist(T1est(T1est>0.01 & T1est < 4.5)*1000, 400); faxis
st(1000*bsxfun(@times, T1est, mask), [0, 2700]); colormap('parula'), colorbar;
title('T1 map (ms)'); faxis

%%
% another display T1
st(1000*bsxfun(@times, T1est, mask), [0, 2700]); colormap('parula'), colorbar;
title('T1 map (ms)'); faxis
stc(bsxfun(@times, proton1, mask)); colormap('parula'), colorbar;
title('mag1 map'); faxis
stc(bsxfun(@times, proton2, mask)); colormap('parula'), colorbar;
title('mag2 map'); faxis


%% 
% T2 fit:
% 300 sec on [256, 256, 6, 32]

[ny, nz, ns, nc, nt] = size(images_TE_gw);
d1 = dimnorm(images_TE_gw, 5);
d2 = d1;
T2mask = d2 > .01*max(d2(:));
echo_times = TE*1e-3;

proton = zeros(ny, nz, ns, nc);
T2est = zeros(ny, nz, ns, nc);

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
        y_cplx =squeeze(images_TE_gw(yy,zz,ss,cc,:));
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
save T2fit.mat
%%
% display T2
T2vals = T2est(T2est~=0);
figure(1); hist(T2vals(T2vals<0.300)*1000, 600); faxis
%
st(1000*bsxfun(@times, T2est, T2mask), [0, 200]); colormap('parula'), colorbar;
title('T2 map(ms)'); faxis
stc(bsxfun(@times, proton, T2mask)); colormap('parula'), colorbar;
title('proton density'); faxis
%%



function kspaces = pfile2kspace(pfile)
    header = GERecon('Pfile.Header',pfile);
    fsePhaseCorrectionBitMask = 131072;
    numRefViews = 0;
    if( bitand(header.RawHeader.data_collect_type1, fsePhaseCorrectionBitMask) )
        numRefViews = header.RawHeader.etl;
    end 
    
    for c = 1:pfile.channels
       kspace = GERecon('Pfile.KSpace',  1,1, c);
            % Discard reference views
       fullYRes = size(kspace, 2);
       yRes = fullYRes - numRefViews;
       kspace = kspace(:, 1:yRes);
       %kspace = kspace(end:-1:1,:);
       kspace(:,2:2:end) = -kspace(:,2:2:end);
       %channelImage = GERecon('Transform', kspace);
       kspaces(:,:,c) = kspace;
       
    end
    
                  
end

                 
