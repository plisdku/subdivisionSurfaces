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

% Special case: a triangle.  It can't be oriented.
% (Any triangles disconnected from the rest of the mesh will also have this
% problem and will end up ignored.  I'm not testing for such a problem.)
if numVertices == 3
    faces = [1 2 3];
    warning('Triangle cannot be oriented.  Picking order [1 2 3].');
end

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

wrap = @(n,N) 1 + mod(n-1,N);

for vv = 1:numVertices
    % putMeFirst([1 7 3 2], 3) returns [3 2 1 7].
    %putMeFirst = @(indexRow, index) circshift(indexRow, ...
    %    [0, 1-find(indexRow == index)]);
    
    if numNeighbors(vv) == 2
        % Because a vertex with two neighbors can't be oriented it's liable
        % to insert a wrongly-oriented face into the output list.  Just let
        % the other vertices that share its (only) triangle deal with the
        % orientation!  Of course if the mesh is merely one triangle then
        % this will fail.
        continue;
    end
    
    for jj = 1:numNeighbors(vv) % for each adjacent face
        ww = VV(vv,jj);
        xx = VV(vv, wrap(jj+1, numNeighbors(vv)));
        
        %if ismember(ww, VV(xx,:)) % if this is a complete triangle
        if find(VV(xx,:) == ww, 1) % if this is a complete triangle
            face = full([vv ww xx]);
            %fprintf('Face is %i %i %i\n', face(1), face(2), face(3));

            if iFace > size(redundantFaces,1)
                warning('Doubling face buffer.');
                redundantFaces(2*iFace,3) = 0; % double its size...
            end
            
            % put the lowest index first.  This will help remove redundant
            % faces later.
            
            if face(1) <= face(2)
                if face(1) <= face(3)
                    redundantFaces(iFace,:) = face;
                else
                    redundantFaces(iFace,:) = face([3 1 2]);
                end
            else
                if face(2) <= face(3)
                    redundantFaces(iFace,:) = face([2 3 1]);
                else
                    redundantFaces(iFace,:) = face([3 1 2]);
                end
            end
            
            %[~,minVertexIndex] = min(face);
            %redundantFaces(iFace,:) = circshift(face, [0, 1-minVertexIndex]);
            %assert(isequal(redundantFaces(iFace,:), putMeFirst(face,min(face))));
            %redundantFaces(iFace,:) = putMeFirst(face, min(face));
            iFace = iFace + 1;
            
        end
    end
end
faces = unique(redundantFaces(1:iFace-1,:), 'rows');


%if size(faces, 1) ~= estimatedNumFaces
%    warning('Expected %i faces from Euler characteristic but got %i faces.', ...
%        estimatedNumFaces, size(faces,1));
%end



