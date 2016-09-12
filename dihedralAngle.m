function beta = dihedralAngle(vv, ww, VV, vertices)
% function beta = dihedralAngle(vv, ww, VV, vertices) 

import VVMesh.*

beta = asin(norm(cross(normalVector(vv, ww, VV, vertices), ...
    normalVector(ww, vv, VV, vertices))));
