function Hv = vertexMeanCurvature(uu, VV, vertices)
% function Hv = vertexMeanCurvature(uu, VV, vertices)
%
% Calculate the integral absolute mean curvature in the area attributed to
% the vertex uu.

import VVMesh.*

valence = nnz(VV(uu,:));

Hv = 0;

for ee = 1:valence
    
    if ~isUnorientedEdgeOnBoundary(uu, VV(uu,ee), VV)
        edgeNorm = norm( vertices(VV(uu,ee),:) - vertices(uu,:) );
        beta = dihedralAngle(uu, VV(uu,ee), VV, vertices);

        Hv = Hv + edgeNorm * abs(beta);
    end
    
end

Hv = 0.25 * Hv;

