% Test!

testEqual = @(a, b) assert(a == b);

% Create a test mesh.  The vertices are arranged in a skew 3x3 square:
%
% 7 8 9
% 4 5 6
% 1 2 3
%
% All vertices are neighbored by their Cartesian neighbors in this array,
% and each square is cut diagonally from upper left to bottom right as
% well.
%

[vertices faces] = flatRegularMesh(3);
VV = fv2vv(faces, vertices);

plotVV(VV, vertices)
view(2);
axis xy image

%% Tests on triangle edges.
% Test both internal and external edges of the mesh (boundaries and
% non-boundaries).  If u-v is an edge of a triangle (there's a triangle to
% the left of u-v), then next(u,v) is the same as prev(u,v).

next = @(u, v) nextInTriangle(u, v, VV);
prev = @(u, v) prevInTriangle(u, v, VV);

testEqual(next(1,2), 4);
testEqual(next(2,3), 5);
testEqual(next(4,1), 2);
testEqual(next(5,2), 3);

testEqual(prev(1,2), 4);
testEqual(prev(2,3), 5);
testEqual(prev(4,1), 2);
testEqual(prev(5,2), 3);

fprintf('Prev and next on triangles, PASSED\n');

%% Tests on boundary edges.
% next(u,v) on a boundary should just be the next edge along that boundary!
% Think of "next" as flopping an edge clockwise around its end.
% Think of "prev" as flopping an edge counterclockwise around its start.

fullSide = [1 2 3 6 9 8 7 4]; % traversing this way: IS NOT a boundary
numBoundaryEdges = numel(fullSide);
wrap = @(n) 1 + mod(n-1, numBoundaryEdges);

for ee = 1:numBoundaryEdges
    v0 = fullSide(ee);
    v1 = fullSide(wrap(ee+1));
    v2 = fullSide(wrap(ee+2));
    
    % Traverse each edge backwards to stay on the outside of the mesh.
    
    testEqual(next(v2,v1), v0);
    testEqual(prev(v1,v0), v2);
end

fprintf('Prev and next on %i boundary edges, PASSED\n', numBoundaryEdges);

%% Distinguishing boundaries
% isEdgeOnBoundary means "is this edge not bounded by a triangle on the
% left."

isEdge = @(u, v) isEdgeOnBoundary(u,v,VV);

for ee = 1:numBoundaryEdges
    v0 = fullSide(ee);
    v1 = fullSide(wrap(ee+1));
    
    testEqual(isEdge(v1,v0), true);
    testEqual(isEdge(v0, v1), false);
end

fprintf('isEdgeOnBoundary around %i outside edges, PASSED\n', ...
    numBoundaryEdges);

%% Interior edges are not boundaries
% Apply isEdgeOnBoundary to some interior edges

v0 = 5; % the center vertex of the space
v1s = VV(v0, 1:nnz(VV(v0,:))); % all neighbors of v0

for v1 = v1s
    testEqual(isEdge(v0, v1), false);
end

fprintf('isEdgeOnBoundary on %i inside edges, PASSED\n', numel(v1s));


