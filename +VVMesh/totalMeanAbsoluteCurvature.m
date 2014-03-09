function Htot = totalMeanAbsoluteCurvature(VV, vertices)
% function Htot = totalMeanAbsoluteCurvature(VV, vertices)

import VVMesh.*

numVertices = size(VV, 1);

valence = numVVNeighbors(VV);

Htot = 0;

for vv = 1:numVertices
for ii = 1:valence(vv)
    ww = full(VV(vv,ii));

    if vv < ww
    if ~isUnorientedEdgeOnBoundary(vv, ww, VV)
        
        edgeLength = norm( vertices(vv,:) - vertices(ww,:) );
        beta = dihedralAngle(vv, ww, VV, vertices);
        
        myH = 0.25*edgeLength*abs(beta);
        
        %fprintf('(%i %i) contributes %f\n', vv, ww, myH);
        
        Htot = Htot + myH;
        
    end
    end

end
end

