function [VV2 vertices2] = extrude(VV, vertices, displacement)
% [VV2 vertices2] = extrude(VV, vertices, displacement)
% 
% Extrude the given structure (presumed to be a mesh with a boundary) along
% the positions in displacement (size [Nd 3]) and triangulate their
% connections.  The final structure will be oriented consistently with the
% original face, if possible.  (Mobius strips need not apply.)
%
% If the mesh has more than one boundary, who knows what will happen!

import VVMesh.*

nearBoundary = VVMesh.boundaryVertices(VV);

N = size(VV,1);
numBoundaryVertices = size(nearBoundary,1);

vertices2 = [vertices; bsxfun(@plus, vertices, displacement)];

% Convert to face-vertex representation.  So much easier that way!
nearFaces = vv2fv(VV);
farFaces = fliplr(nearFaces + N);

% Link up boundary vertices.

wrap = @(n) 1 + mod(n-1, N);

sideFaces = zeros(numBoundaryVertices*2, 3);
sideFaces = [nearBoundary, wrap(nearBoundary-1), wrap(nearBoundary-1)+N; ...
    wrap(nearBoundary-1)+N, nearBoundary+N, nearBoundary];

faces = [nearFaces; sideFaces; farFaces];

%checkFaces(faces)

VV2 = fv2vv(faces, vertices2);

if ~isOrientedConsistently(VV2)
    error(['Orientation of extruded mesh is inconsistent.  ' ...
        'Did you extrude upside-down?']);
end

% Fix the orientation to agree with the original face.

