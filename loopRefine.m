function [VV2 T vertices2 inheritance] = loopRefine(VV, vertices)
% [outVV T refinedVertices inheritance] = loopRefine(VV, vertices)
%
% Mesh refinement step in Loop subdivision.  Does not perturb
% vertices--only handles the refinement step!  (New mid-edge vertices are
% just at the midpoints of edges, then.)
%
% outVV is the refined vertex-vertex topology.
%
% T is the matrix of weights to perturb the vertices by.  To accomplish
% the perturbation you can set vertices2 = T*vertices yourself.
%
% vertices2 contains an unsmoothed refined set of vertices.
%
% inheritance is a sparse matrix containing the pairs of original vertex
% indices that gave rise to each split vertex.  inheritance(vNew) = [v1 v2]
% means that vNew split the edge from v1 to v2.
% 

numOriginalVertices = size(VV,1);
assert(size(vertices,1) == numOriginalVertices);

%if nargin < 3
%    splitAdjacentEdgeFlags = ones(numOriginalVertices, 1);
%end

numNeighbors = numVVNeighbors(VV);
totalNeighbors = sum(numNeighbors);
assert(rem(totalNeighbors,2) == 0);
numEdges = totalNeighbors/2;

%% Helper functions.

% Helper methods for incrementing indices in wraparound fashion, 1-based.
circPrev = @(n, nMax) 1 + mod(n-2, nMax);
circNext = @(n, nMax) 1 + mod(n, nMax); % (n+1) - 1 == n.
findIndex = @(iFindMe, indices) find(indices == iFindMe, 1);

next = @(v0, v1, myVV) myVV(v1, circPrev(...
    findIndex(v0, myVV(v1,:)), nnz(myVV(v1,:))));
prev = @(v0, v1, myVV) myVV(v0, circNext(...
    findIndex(v1, myVV(v0,:)), nnz(myVV(v0,:))));
isBoundary = @(u,v) next(u,v,VV) ~= prev(u,v,VV);

vertexOnBoundary = @(v) any(arrayfun(@(w) full(isBoundary(v,w)), ...
    full(VV(v, 1:numNeighbors(v)))));

%% Add a vertex at the midpoint of each edge!

vertices2 = vertices;
VV2 = VV;

inheritance = spalloc(numOriginalVertices + numEdges, 2, ...
    2*numEdges);

vNew = numOriginalVertices + 1; % index of next new vertex!

% Iterate over edges:
%iSplit = find(splitAdjacentEdgeFlags);
for v0 = 1:numOriginalVertices
for jj = 1:numNeighbors(v0)
    v1 = VV(v0,jj);
        
    if v0 < v1 % this step prevents processing an edge twice
    %if all(splitAdjacentEdgeFlags([v0 v1]))
        
        % Split the edge!
        vertices2(vNew,:) = 0.5*(vertices(v0,:) + vertices(v1,:));
        
        % Edge is now (u w v).
        %  u replaces v with w
        %  v replaces u with w
        %  w adds u and v (either order is ok, but I'll go v1-v0)
        
        VV2(v0, VV2(v0,:) == v1) = vNew; % v0 replaces v1 with vMid
        VV2(v1, VV2(v1,:) == v0) = vNew; % v1 replaces v0 with vMid
        VV2(vNew, 1) = v1; % vMid adds v0
        VV2(vNew, 2) = v0; % vMid adds v1
        
        % and that concludes the addition of vertices; it's no longer a
        % triangular mesh though!
        
        inheritance(vNew, [1 2]) = [v0 v1];
        
        vNew = vNew + 1;
    %end
    end
end
end

numNeighbors2 = numNeighbors;
numNeighbors2(numOriginalVertices+1:size(vertices2,1)) = 2;

checkSymmetry(VV2);

%% Now I can allocate the weight matrix.
% The number of nonnzeros should be
%  numOriginalVertices + totalNeighbors  ... (for the old vertices)
% ... + 4*numNewVertices (for the new vertices)
%
% but it might end up smaller if the mesh has a boundary.

% First handle the original vertices.

numVertices = size(vertices2, 1);
numNewVertices = numVertices - numOriginalVertices;

T = spalloc(numVertices, numOriginalVertices, ...
    numOriginalVertices + totalNeighbors + 4*numNewVertices);

for vOld = 1:numOriginalVertices
    vNeighbors = VV(vOld, VV(vOld,:) ~= 0);
    
    neighborWeight = loopWeights(numel(vNeighbors));
    selfWeight = 1 - neighborWeight*numel(vNeighbors);
    
    if vertexOnBoundary(vOld)
        T(vOld, vOld) = 1.0; % don't move off boundary.
    else
        T(vOld, vOld) = selfWeight;
        T(vOld, vNeighbors) = neighborWeight;
    end
end

%% Now add connections to VV2 to make it a triangulation again.


printVert = @(ii) fprintf('Vertex %i is at %i %i %i\n', full(ii), ...
    full(vertices2(ii,1)), full(vertices2(ii,2)), full(vertices2(ii,3)));

isOld = @(vIndex) vIndex <= numOriginalVertices;

% For each NEW vertex:
for vNew = numOriginalVertices+1:size(vertices2,1)
    
    %printVert(ii);
    
    boundaryFlag = false;
    
    v1 = VV2(vNew,1);
    v0 = VV2(vNew,2);
    
%     fprintf('From ');
%     printVert(v0);
%     fprintf('Through ');
%     printVert(vNew);
%     fprintf('To ');
%     printVert(v1);
    
    % First handle the triangle that was bounded by (v0 v1)
    if isBoundary(v0, v1)
        neighbors01 = v1;
        boundaryFlag = true;
        
%         fprintf('%i to %i was an outer boundary.\n', full(v0), full(v1));
%         fprintf('v0: \n');
%         disp(full(VV2(v0, 1:4)));
%         fprintf('vNew: \n');
%         disp(full(VV2(vNew, 1:4)));
%         fprintf('v1: \n');
%         disp(full(VV2(v1, 1:4)));
        % v0-vNew-v1 is a boundary, meaning it's empty to the left.
        % do nothing.
    else
        vMidNext = next(vNew, v1, VV2);
        vMidPrev = prev(v0, vNew, VV2);
        
        if isOld(vMidNext) % it's a T-vertex
            neighbors01 = [v1 vMidNext];
            if vMidNext ~= vMidPrev
                error('T vertex has extra neighbor!');
            end
            
            N = numNeighbors2(vMidNext);
            here = findIndex(v1, VV2(vMidNext,:));
            VV2(vMidNext, 1:N+1) = [VV2(vMidNext, 1:here-1), vNew, ...
                VV2(vMidNext, here:N)];
            numNeighbors2(vMidNext) = N + 1;
        else
            neighbors01 = [v1 vMidNext vMidPrev];
        end
        
    end
    
    % Next handle the triangle that was bounded by (v1 v0)
    if isBoundary(v1, v0)
        neighbors10 = v0;
        boundaryFlag = true;
        
%         fprintf('%i to %i was an outer boundary.\n', full(v1), full(v0));
%         fprintf('v1: \n');
%         disp(full(VV2(v1, 1:4)));
%         fprintf('vNew: \n');
%         disp(full(VV2(vNew, 1:4)));
%         fprintf('v0: \n');
%         disp(full(VV2(v0, 1:4)));
        % v0-vNew-v1 is a boundary, meaning it's empty to the left.
        % do nothing.
    else
        
        vMidNext = next(vNew, v0, VV2);
        vMidPrev = prev(v1, vNew, VV2);
        
        if isOld(vMidNext) % it's a T-vertex
            neighbors10 = [v0 vMidNext];
            if vMidNext ~= vMidPrev
                error('T vertex has extra neighbor!');
            end
            
            N = numNeighbors2(vMidNext);
            here = findIndex(v0, VV2(vMidNext,:));
            VV2(vMidNext, 1:N+1) = [VV2(vMidNext, 1:here-1), vNew, ...
                VV2(vMidNext, here:N)];
            numNeighbors2(vMidNext) = N + 1;
        else
            neighbors10 = [v0 vMidNext vMidPrev];
        end
        
    end
    
    assert(~ismember(vNew, neighbors01));
    assert(~ismember(vNew, neighbors10));
    assert(numel(unique(neighbors01)) == numel(neighbors01));
    assert(numel(unique(neighbors10)) == numel(neighbors10));
    
    % All neighbors for the vertex!!
    VV2(vNew, 1:(numel(neighbors01) + numel(neighbors10))) = ...
        [neighbors01 neighbors10];
    
    if numel(unique([neighbors01 neighbors10])) <= 2
        error('Fuck.');
    end
    
    % Now the T-matrix:
    if boundaryFlag
        T(vNew, [v0 v1]) = 0.5;
    else
        T(vNew, [v0 v1]) = 0.375;
        T(vNew, [next(v1, v0, VV) next(v0, v1, VV)]) = 0.125;
    end
end

%% Tests!

getVV2 = @(v) full(VV2(v, 1:nnz(VV2(v,:))));
checkSymmetry(VV2);

Tsums = sum(T,2);

if ~all(abs(Tsums - 1.0) < 1e-6)
    error('T seems wrong.\n');
end


