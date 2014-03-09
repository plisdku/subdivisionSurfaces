function checkSymmetry(VV)

import VVMesh.*

A = vv2adjacency(VV);

if ~isequal(A, A')
    error('Not symmetrical.');
end