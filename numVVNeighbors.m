function nn = numVVNeighbors(VV)
% nn = numNeighbors(VV)
% Calculate number of neighbor vertices for each vertex in a VV mesh

import VVMesh.*

numVertices = size(VV,1);
nn = zeros(numVertices,1);

for ii = 1:numVertices
    nn(ii) = nnz(VV(ii,:));
end
