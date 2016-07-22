function [range, tot] = st(in, range, shape, nip, sp, colors, figname, titlename)
% Inputs:
%     - in: matrix or cell array
%           if cell array, elements must have the same first and second
%           dimensions
%     - range: window range, can be one of the following
%             # []: default = [min(:) max(:)]
%             # [i 1]: [min(:) 1]
%             # [1 i]: [1 max(:)]
%     - shape: one of the following
%                + []: default (aspect ratio = 1)
%                + array specifying dimension of montage. EX:
%                       # [4 2]
%                       # [2 i] first dim is 2, second dim is calc.
%                + desired aspect ratio = x/y
%     - nip: one of the following
%                + 0: do nothing
%                + 1: show image
%                + 2: print to file
%                + []: default (1)
%     - colors: [fillc spc], spc is color for spacers, fillc is color for fillers
%                + []: default ([0 1], i.e. white and black)
%     - sp: [vsp hsp], vsp is vertical spacing, hsp is horizontal spacing
%               + []: default ([5 5])
%     - figname: name for figure
%                + []: default (no figure)
% Outputs:
%     - range: window range
%     - tot: tiled image
% Notes:
%    - wrg.m might be useful
% D. Bahri

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
    in = cellfun(@(x) reshape(x, [size(x,1) size(x,2) numel(x)/(size(x,1)*size(x,2))]), in, 'UniformOutput', false);
    in = double(cat(3, in{:}));
else
    in = double(in);
end

if isempty(range)
    range = [min(in(:)) max(in(:))];
elseif ~isreal(range(1))
    range = [min(in(:)) range(2)];
elseif ~isreal(range(2))
    range = [range(1) max(in(:))];
end

spval = colors(2)*(range(2)-range(1)) + range(1);
fillval = colors(1)*(range(2)-range(1)) + range(1);
vsp = sp(1);
hsp = sp(2);

nx = size(in,1);
ny = size(in,2);
nz = numel(in)/(nx*ny);
in = reshape(in, [nx ny nz]);

if length(shape) == 1
    % fit to aspect ratio as best as possible
    y = ceil(sqrt(nz/shape));
    x = ceil(y*shape);
elseif ~isreal(shape(1))
    y = shape(2);
    x = ceil(nz/y);
elseif ~isreal(shape(2))
    x = shape(1);
    y = ceil(nz/x);
else
    x = shape(1);
    y = shape(2);
end

tot = [];
for j = 0:y-1
    tmp = [];
    for i = 0:x-1
        idx = i+j*x+1;
        if idx > nz
            mat = fillval*ones([nx ny]);
        else
            mat = in(:,:,idx);
        end
        if i ~= x-1
            mat = vertcat(mat,spval*ones([vsp ny]));
            if j ~= y-1
                mat = horzcat(mat,spval*ones([nx+vsp hsp]));
            end
        elseif j ~= y-1
            mat = horzcat(mat,spval*ones([nx hsp]));
        end
        tmp = vertcat(tmp,mat);
    end
    tot = horzcat(tot,tmp);
end

switch nip
    case 0
        return;
    case 1
        if isempty(figname)
            figure; imshow(tot,range);
            impixelinfo;
            title(titlename);
            return;
        else
            figure('Name',figname); imshow(tot,range);
            impixelinfo;
            title(titlename);
            return;
        end
    case 2
        if isempty(figname)
            error('must have figname');
        else
            %export_fig(figname,'-pdf','-q101');
            h = figure; imshow(tot,range);
            title(titlename);
            pdfcrop(h);
            saveas(h,figname,'pdf');
%             close(h)
            % needs texlive
%             unix(['/usr/texbin/pdfcrop ' figname '.pdf ' figname '.pdf']);
        end
end



