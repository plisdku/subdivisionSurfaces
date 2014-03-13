function yesNo = isOrientedConsistently(VV)
% yesNo = isOrientedConsistently(VV)
%
% Check that all vertices in VV are oriented consistently with their
% neighbors.

import VVMesh.*

for vv = 1:size(VV,1)
    for ww = nonzeros(VV(vv,:))'
        if ~areSameOrientation(vv, ww, VV)
            yesNo = false;
            return
        end
    end
end

yesNo = true;