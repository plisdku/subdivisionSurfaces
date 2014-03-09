function faces = vv2fv(VV, numNeighbors)
% faces = vv2fv(VV)
% faces = vv2fv(VV, numNeighbors)
%
% Return the face index table for a mesh in vertex-vertex representation.
%
% VV = [Nv maxNumNeighbors] sparse array of vertex indices.
% 
% numNeighbors = [Nv 1] array of how many neighboring vertices each vertex
% has.  Optional (can be inferred from VV).
%
% faces = [Nf 3] array of vertex indices.
%
% maxNumNeighbors is at least the maximum number of neighbors of any vertex
% in the mesh, but could be larger.
%
%

import VVMesh.*

numVertices = size(VV,1);

if nargin == 1
    numNeighbors = zeros(numVertices, 1);
    
    for ii = 1:numVertices
        numNeighbors(ii) = nnz(VV(ii,:));
    end    
end

% Use Euler characteristic to determine num faces: V - E + F = 2 for
% polyhedra.
totalNeighbors = sum(numNeighbors);
assert(rem(totalNeighbors,2) == 0);
numEdges = totalNeighbors/2;

estimatedNumFaces = 2 - numVertices + numEdges;

redundantFaces = zeros(3*estimatedNumFaces, 3);

iFace = 1;

for ii = 1:numVertices
    % putMeFirst([1 7 3 2], 3) returns [3 2 1 7].
    putMeFirst = @(indexRow, index) circshift(indexRow, ...
        [0, 1-find(indexRow == index)]);
    
    for jj = 1:numNeighbors(ii)-1 % for each adjacent face
        face = [ii, VV(ii,jj), VV(ii,jj+1)];
        
        % put the lowest index first.  This will help remove redundant
        % faces later.
        redundantFaces(iFace,:) = putMeFirst(face, min(face));
        iFace = iFace + 1;
    end
end

faces = unique(redundantFaces(1:iFace-1,:), 'rows');


%if size(faces, 1) ~= estimatedNumFaces
%    warning('Expected %i faces from Euler characteristic but got %i faces.', ...
%        estimatedNumFaces, size(faces,1));
%end



