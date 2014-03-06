%function VV2 = optimizeVVFaces(VV, vertices)
% VV2 = optimizeVVFaces(VV, vertices)
%
% Iteratively flip diagonals of split quads to make triangles more uniform.

numVertices = size(VV,1);

valence = numVVNeighbors(VV);
numEdges = sum(valence)/2;

vertexOnBoundary = @(v) any(arrayfun(@(w) full(isEdgeOnBoundary(v,w,VV)), ...
    full(VV(v, 1:valence(v)))));

%%

edges = zeros(numEdges, 2);
diags = zeros(numEdges, 2);
scores = zeros(numEdges, 1);

currEdge = 1;

for vv = 1:numVertices
    
    for ii = 1:valence(vv)
    ww = VV(vv,ii);
    if vv < ww && ~isUnorientedEdgeOnBoundary(vv,ww,VV)
        
        % Get the other two vertices of the split quad
        xx = nextInTriangle(vv,ww,VV);
        yy = nextInTriangle(ww,vv,VV);
        
        % Silly idea: rate by the ratio of diagonal lengths.
        
        score = norm(vertices(ww,:) - vertices(vv,:)) / ...
            norm(vertices(xx,:) - vertices(yy,:));
        
        edges(currEdge,:) = [vv ww];
        diags(currEdge,:) = [xx yy];
        scores(currEdge,:) = score;
        
        currEdge = currEdge + 1;
    end
    end
    
end

%%