function labelVertices(VV, vertices, sel, varargin)
% labelVertices(VV, vertices, sel)
%
% Write vertex indices on a 3D plot

numVertices = size(vertices,1);

if nargin < 3
    sel = [];
end

if isempty(sel)
    sel = 1:numVertices;
elseif islogical(sel)
    assert(numel(sel) == numVertices);
    sel = find(sel);
end

for vv = reshape(sel, 1, [])
    text(vertices(vv,1), vertices(vv,2), vertices(vv,3), ...
        ['  ', num2str(vv)], varargin{:});
end
