function [faceEdgesUnique, faceEdgeOrientations, edgeVerticesUnique] = getEdges(faces)

numFaces = size(faces, 1);

% First make a raw edge list
faceEdges = reshape(1:3*numFaces, numFaces, 3);

% Now make a list of vertices for each edge
edgeVertices = [faces(:, 1:2); faces(:,2:3); faces(:,[3,1])];

% Now make the unique ones
edgeVerticesSorted = sort(edgeVertices, 2);
[cc, iaa, icc] = unique(edgeVerticesSorted, 'rows');
%%
edgeVerticesUnique = edgeVertices(iaa,:);
faceEdgesUnique = reshape(icc(faceEdges), numFaces, 3);

%% Get orientation of each edge-use of each face!
faceEdgeOrientations = 0*faceEdgesUnique;

for iFace = 1:numFaces
    face = faces(iFace,:);
    for iEdge = 1:3
        e0 = edgeVerticesUnique(faceEdgesUnique(iFace,iEdge),:);
        v0 = faces(iFace,iEdge);

        if (e0(1) == v0)
            faceEdgeOrientations(iFace,iEdge) = 1;
        else
            assert(e0(2) == v0)
            faceEdgeOrientations(iFace,iEdge) = -1;
        end
    end
end

