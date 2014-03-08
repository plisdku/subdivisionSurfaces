function nv = normalVector(vv, ww, VV, vertices)
% function nv = normalVector(vv, ww, VV, vertices)

if isEdgeOnBoundary(vv, ww, VV)
    error('Boundary edge does not refer to triangle at all');
end

xx = nextInTriangle(vv, ww, VV);

up = cross(vertices(ww,:) - vertices(vv,:), ...
    vertices(xx,:) - vertices(vv,:));

nv = up / norm(up);

