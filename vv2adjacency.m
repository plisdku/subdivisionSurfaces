function A = vv2adjacency(VV)
% A = vv2adjacency
% 
% Create vertex-vertex adjacency matrix.

import VVMesh.*

[rr cc vv] = find(VV);

numVertices = size(VV,1);
A = sparse(rr, vv, ones(size(vv)), numVertices, numVertices);