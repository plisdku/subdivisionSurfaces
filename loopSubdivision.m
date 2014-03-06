function [VV2 vertices2] = loopSubdivision(VV, vertices)

numOriginalVertices = size(vertices,1);

%if nargin < 3
%    splitAdjacentEdgeFlags = ones(numOriginalVertices, 1);
%end

numNeighbors = numVVNeighbors(VV);

[VV2, T] = loopRefine(VV, vertices);
% ok that was the easy part.

vertices2 = T * vertices;
