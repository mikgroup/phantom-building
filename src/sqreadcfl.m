function data = sqreadcfl(filenameBase)
% function data = sqreadcfl(filenameBase)
%
% Read in recon data stored in filenameBase.cfl (complex float)
% based on dimensions stored in filenameBase.hdr.
% JT
% squeeze degenerate dimensions
[a, b, ~] = fileparts(filenameBase);
data = squeeze(readcfl([a, '/', b]));
end
