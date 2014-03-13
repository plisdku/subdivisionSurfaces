function VV = polygon(points)
% VV = polygon(points)
% 
% Triangulate the interior of the region bounded by points.
% points may be 2d or 3d.

import VVMesh.*

nDims = size(points,2);
numVertices = size(points,1);

if nDims == 2
    constraints = [(1:numVertices)', ([2:numVertices,1]')];
    dt = DelaunayTri(points(:,1), points(:,2), constraints);
    
    in = inOutStatus(dt);
    faces = dt.Triangulation(in,:);
    
else
    % 3D points.  I'll assume they're supposed to be essentially planar and
    % flatten out the extra dimension in a sensible way.
    
    center = mean(points);
    ptsCentered = bsxfun(@plus, points, -center);
    
    [~,~,V] = svd(ptsCentered'*ptsCentered);
    pts2d = ptsCentered * V(:,1:2);
    
    constraints = [(1:numVertices)', ([2:numVertices,1]')];
    dt = DelaunayTri(pts2d(:,1), pts2d(:,2), constraints);
    
    in = inOutStatus(dt);
    faces = dt.Triangulation(in,:);
end


VV = fv2vv(faces);
