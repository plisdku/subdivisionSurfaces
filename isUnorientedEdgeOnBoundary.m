function onBoundary = isUnorientedEdgeOnBoundary(u, v, VV)
% onBoundary = isUnorientedEdgeOnBoundary(v0, v1, VV)
%
% Determine whether u-v is bounded on both sides by triangles.

import VVMesh.*

onBoundary = (nextInTriangle(u,v,VV) ~= prevInTriangle(u,v,VV)) || ...
    (nextInTriangle(v,u,VV) ~= prevInTriangle(v,u,VV));