function vNext = nextInTriangle(v0, v1, VV)
% vNext = nextInTriangle(v0, v1, VV)
%
% If v0-v1 is an edge of a triangle, find the next vertex v2 that completes
% the triangle, in counterclockwise orientation.
%
% This is equivalent to traversing from v0 to v1 and then choosing the
% leftmost edge to follow next.

import VVMesh.*

circPrev = @(n, nMax) 1 + mod(n-2, nMax);
findIndex = @(iFindMe, indices) find(indices == iFindMe, 1);

vNext = full(VV(v1, circPrev(findIndex(v0, VV(v1,:)), nnz(VV(v1,:)))));
