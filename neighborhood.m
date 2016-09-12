function ww = neighborhood(vv, VV, radius, minCount)
% ww = neighborhood(vv, VV)
%
% Return all vertices in the 1-neighborhood of vertices vv.
%
% neighborhood(vv, VV, r) gives the r-neighborhood of vertices vv.  It's
% done by recursion so it's not going to be that fast.
%
% neighborhood(vv, VV, r, minCount) returns vertices in the r-neighborhood
% of vertices vv which neighbor at least minCount of the vertices in vv.
% The default behavior is as though minCount == 1.

import VVMesh.*

if nargin < 3
    radius = 1;
else
    assert(radius == round(radius)); % radius should be an integer
    assert(radius >= 0);
end

if nargin < 4
    minCount = 1;
end

numVertices = size(VV,1);

if radius == 1
    ww = find(histc(nonzeros(VV(vv,:)),1:numVertices) >= minCount);
elseif radius > 1
    ww = find(histc(nonzeros(VV(neighborhood(vv, VV, radius-1),:)),...
        1:numVertices) >= minCount);
elseif radius == 0
    ww = full(vv);
else
    error('This should never happen.');
end
