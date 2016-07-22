function [tot, range] = stc(in, range, shape, nip, sp, colors, figname, titlename)
% Overloaded for abs of st. See st.m for usage
% JT

% defaults
if nargin < 8
    titlename = '';
end
if nargin < 7
    figname = [];
end

if nargin < 6 || isempty(colors)
    colors = [0 1];
end

if nargin < 5 || isempty(sp)
    sp = [5 5];
end

if nargin < 4 || isempty(nip)
    nip = 1;
end

if nargin < 3 || isempty(shape)
    shape = 1;
end

if nargin < 2
    range = [];
end

if iscell(in)
    for i=1:length(in)
        in{i} = abs(in{i});
    end
else
    in = abs(in);
end

[tot, range] = st(in, range, shape, nip, sp, colors, figname, titlename);
end


