function VV = polygon(points, hole)
% VV = polygon(points)
% 
% Triangulate the interior of the region bounded by points.
% points may be 2d or 3d; if 3D the direction of least variance will be
% disregarded.
%
% VV = polygon(points, holePoints)
%
% Triangulate the interior of a polygon that has a hole.

    import VVMesh.*

    nDims = size(points,2);
    
    if nargin == 2
        nOuter = size(points,1);
        nHole = size(hole, 1);
        nVertices = nOuter + nHole;
        
        outerConstraint = [ (1:nOuter)', [2:nOuter,1]' ];
        innerConstraint = nOuter + [ (1:nHole)', [2:nHole,1]' ];
        
        constraints = [outerConstraint; innerConstraint];

        allPoints = [points; hole];
    else
        nVertices = size(points,1);
        constraints = [(1:nVertices)', ([2:nVertices,1]')];

        allPoints = points;
    end

    if nDims == 3
        center = mean(allPoints);
        ptsCentered = bsxfun(@plus, allPoints, -center);

        [~,~,V] = svd(ptsCentered'*ptsCentered);
        pts2d = ptsCentered * V(:,1:2);
        
        faces = myTris(pts2d, constraints);
    else
        faces = myTris(allPoints, constraints);
    end
    
    VV = fv2vv(faces);
end




function faces = myTris(pts2d, boundaryEdges)
    dt = DelaunayTri(pts2d(:,1), pts2d(:,2), boundaryEdges);
    in = inOutStatus(dt);
    faces = dt.Triangulation(in,:);
end

