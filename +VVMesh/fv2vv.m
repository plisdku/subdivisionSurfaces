function VV = fv2vv(faces, vertices)
% vv = fv2vv(faces)
% vv = fv2vv(faces, vertices)
%
% Return the vertex neighbor table vv for a mesh in face-vertex
% representation.
%
% faces = [Nf 3] array of vertex indices
% vv = [Nv maxNumNeighbors] sparse array of vertex indices.
%
% maxNumNeighbors is at least the maximum number of neighbors of any vertex
% in the mesh, but could be larger.
%

import VVMesh.*

if nargin == 1
    numVertices = numel(unique(faces(:)));
else
    numVertices = size(vertices, 1);
end

numFaces = size(faces, 1);

assert(max(faces(:)) <= numVertices);
assert(min(faces(:)) >= 1);

% Euler characteristic: V - E + F = 2 for polyhedra.
% Use to guess number of edges.
% This tells us how to guess the number of VV neighbors.
% Should help with the sparse matrix writing.

estimatedNumEdges = numVertices + numFaces - 2;
estimatedTotalNeighbors = 2*estimatedNumEdges; % two vertices share an edge

MANY_VERTEX_NEIGHBORS = 100;

VV = spalloc(numVertices, MANY_VERTEX_NEIGHBORS, estimatedTotalNeighbors);

%% Now starts the work of the awesome vertex man.  Yeah!

% These sparse matrix buffers store edges that are opposite to each vertex
% in various faces.  So, vertex ii is opposite edges from vNext(ii,:) to
% vPrev(ii,:).  Use nVertOppositeEdges to record the total number of
% opposite edges.
vPrev = VV;
vNext = VV;
nVertOppositeEdges = zeros(numVertices, 1);

% Number of neighbor vertices recorded in the data structure.
% Used for fast, correct indexing into the sparse VV matrix.
numVertNeighbors = zeros(numVertices,1);

% Update all v1 with neighbors:

for ff = 1:numFaces
    
    face = faces(ff,:);
    
    % Update v1 with neighbors:
    nVertOppositeEdges(face(1)) = nVertOppositeEdges(face(1)) + 1;
    vNext(face(1), nVertOppositeEdges(face(1))) = face(2);
    vPrev(face(1), nVertOppositeEdges(face(1))) = face(3);
    
    % Update v2 with neighbors:
    nVertOppositeEdges(face(2)) = nVertOppositeEdges(face(2)) + 1;
    vNext(face(2), nVertOppositeEdges(face(2))) = face(3);
    vPrev(face(2), nVertOppositeEdges(face(2))) = face(1);
    
    % Update v3 with neighbors:
    nVertOppositeEdges(face(3)) = nVertOppositeEdges(face(3)) + 1;
    vNext(face(3), nVertOppositeEdges(face(3))) = face(1);
    vPrev(face(3), nVertOppositeEdges(face(3))) = face(2);
    
end

%% Chain them together somehow.

% For each vertex, order its neighbors counterclockwise.
% Each edge in the ring going counterclockwise around the vertex is
% a segment from a vNext to a vPrev.  When one of these edges ends on the
% beginning of another edge then I have an ordering of a few of those
% vertices.  So I need to find that ordering!

% There is a face in which vNext(vCurr, someFace) is the next vert, and
% vPrev(vCurr, someFace) is the previous vert.  That's why this works.

for iVert = 1:numVertices
    
    %fprintf('iVert = %i\n', iVert);
    numLoopEdges = nVertOppositeEdges(iVert);
    
    edges = [vNext(iVert,1:numLoopEdges); vPrev(iVert,1:numLoopEdges)];
    % edges is [2 N_edges] in size.  Now perform the ordering.
    
    % Here's the version for surfaces that may have boundaries.
    % We need to find the "first" edge, if there is one.
    % Look for a vertex that is present in edges(1,:) but not in
    % edges(2,:).
    
    [~, startingEdge] = setdiff(edges(1,:), edges(2,:));
    if isempty(startingEdge) % vertex is not on a boundary
        startingEdge = 1;
        numNeighborVerts = numLoopEdges;
    else
        numNeighborVerts = numLoopEdges + numel(startingEdge);
    end
    
    if numel(startingEdge) > 1
        error('Presently not considering meshes where boundaries touch.');
    end
    
    vNeighbors = zeros(1, numNeighborVerts);
    vNeighbors(1:2) = edges(1:2, startingEdge);
    
    for nn = 3:numNeighborVerts
        currentEdge = find(edges(1,:) == vNeighbors(nn-1));
        
        assert(numel(currentEdge) == 1);
        
        vNeighbors(nn) = edges(2,currentEdge);
    end
    
    VV(iVert,1:numNeighborVerts) = vNeighbors;
end






