function v0 = prevInTriangle(v1, v2, VV)
% vPrev = prevInTriangle(v0, v1, VV)
%
% If v1-v2 is an edge of a triangle, find the previous vertex v2 that
% completes the triangle, in counterclockwise orientation.
%
% This is equivalent to traversing from v2 to v1 and then choosing the
% rightmost edge to follow next.

circNext = @(n, nMax) 1 + mod(n, nMax); % (n+1) - 1 == n.
findIndex = @(iFindMe, indices) find(indices == iFindMe, 1);

v0 = full(VV(v1, circNext(findIndex(v2, VV(v1,:)), nnz(VV(v1,:)))));
