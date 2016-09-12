function onBoundary = isEdgeOnBoundary(u, v, VV)
% onBoundary = isEdgeOnBoundary(u, v, VV)
%
% Determine whether the edge u-v is not bounded on the left by a triangle.
%
% This test will work for all boundaries except triangular boundaries,
% since the inside and outside of a triangle both look like a triangle.

import VVMesh.*

onBoundary = nextInTriangle(u,v,VV) ~= prevInTriangle(u,v,VV);
