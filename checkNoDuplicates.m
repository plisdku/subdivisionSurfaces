function checkNoDuplicates(VV)
% checkNoDuplicates(VV) will fail if any vertices have duplicate neighbors.

import VVMesh.*

numVertices = size(VV,1);

for vv = 1:numVertices
    valence = nnz(VV(vv,:));
    numUniques = numel(unique(VV(vv,:))) - 1; % ignore 0, so subtract 1
    
    if valence ~= numUniques
        error('Got duplicate vertex neighbors.');
    end
end