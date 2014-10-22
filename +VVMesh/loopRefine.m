function [VV2, vertices2, T, perturbFlags, T_refine] = loopRefine(VV, vertices, ...
    refineVertices)
% [outVV outVertices T perturbFlags] = loopRefine(VV, vertices)
%
% Mesh refinement step in Loop subdivision.  Does not perturb
% vertices--only handles the refinement step!  (New mid-edge vertices are
% just at the midpoints of edges, then.)
%
% outVV is the refined vertex-vertex topology.
%
% outVertices contains an unsmoothed refined set of vertices.
%
% T is the matrix of weights to perturb the vertices by.  To accomplish
% the perturbation you can set vertices2 = T*vertices yourself.
%
% loopRefine(VV, vertices, refineVertices) will restrict the refinement
% region, as for adaptive refinement.
%
% perturbFlags will be valued true for each vertex that should be allowed
% to move in subdivision.  This is for adaptive refinement with creases.
% 

import VVMesh.*

checkForErrors = false;

numOriginalVertices = size(VV,1);
assert(size(vertices,1) == numOriginalVertices);

if nargin < 3
    refineVertices = 1:numOriginalVertices;
    refineNeighborhood = refineVertices;
else
    refineNeighborhood = union(refineVertices, ...
        neighborhood(refineVertices, VV));
end

splitAdjacentEdgeFlags = zeros(numOriginalVertices,1);
splitAdjacentEdgeFlags(refineNeighborhood) = true;

numNeighbors = numVVNeighbors(VV);
totalNeighbors = sum(numNeighbors);
assert(rem(totalNeighbors,2) == 0);
numEdges = totalNeighbors/2;


%% Add a vertex at the midpoint of each edge!

[vertices2, VV2] = splitEdges(VV, vertices, splitAdjacentEdgeFlags, numNeighbors);

numNeighbors2 = numNeighbors;
numNeighbors2(numOriginalVertices+1:size(vertices2,1)) = 2;

checkSymmetry(VV2);

%% Now I can allocate the weight matrix.
% The number of nonnzeros should be
%  numOriginalVertices + totalNeighbors  ... (for the old vertices)
% ... + 4*numNewVertices (for the new vertices)
%
% but it might end up smaller if the mesh has a boundary or when it's only
% partially refined.

% First handle the original vertices.

numVertices = size(vertices2, 1);
numNewVertices = numVertices - numOriginalVertices;
T = spalloc(numVertices, numOriginalVertices, ...
    numOriginalVertices + totalNeighbors + 4*numNewVertices);

% The matrix that achieves mesh refinement is look like this!
T_refine = sparse(1:numOriginalVertices, 1:numOriginalVertices, ...
    ones(numOriginalVertices,1), numVertices, numOriginalVertices);

perturbedVertices = union(refineVertices, ...
    neighborhood(refineVertices, VV2, 1, 2));
perturbFlags = zeros(numVertices,1);
perturbFlags(perturbedVertices) = true;

% A helper function.
vertexOnBoundary = @(v) any(arrayfun(@(w) full(isEdgeOnBoundary(v,w,VV)), ...
    full(VV(v, 1:numNeighbors(v)))));

for vOld = 1:numOriginalVertices
    if vertexOnBoundary(vOld) || ~perturbFlags(vOld)
        T(vOld, vOld) = 1.0; % don't move off boundary.
    else
        vNeighbors = VV(vOld, VV(vOld,:) ~= 0);

        neighborWeight = loopWeights(numel(vNeighbors));
        selfWeight = 1 - neighborWeight*numel(vNeighbors);

        T(vOld, vOld) = selfWeight;
        T(vOld, vNeighbors) = neighborWeight;
    end
end

%% Now add connections to VV2 to make it a triangulation again.

printVert = @(ii) fprintf('Vertex %i is at %i %i %i\n', full(ii), ...
    full(vertices2(ii,1)), full(vertices2(ii,2)), full(vertices2(ii,3)));

isOld = @(vIndex) vIndex <= numOriginalVertices;

% If I were to use sparse() to speed up some of this stuff, which would be
% a big improvement, I could start like this... it seems like a pain...
%numNewVertices = size(vertices2,1) - numOriginalVertices;
%newRows = zeros(numNewVertices*4,1);
%newCols = zeros(numNewVertices*4,1);
%newVals = zeros(numNewVertices*4,1);
%iNew = 1;

% For each NEW vertex:
for vNew = numOriginalVertices+1:size(vertices2,1)
    
    %if mod(vNew, 1000) == 0
    %    fprintf('vNew = %i\n', vNew);
    %end
    
    %if vNew > size(vertices2,1) || size(vertices2,1) > 17000
    %    keyboard
    %end
    
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
    if isEdgeOnBoundary(v0, v1, VV)
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
        vMidNext = nextInTriangle(vNew, v1, VV2);
        vMidPrev = prevInTriangle(v0, vNew, VV2);
        
        if isOld(vMidNext) % it's a T-vertex
            neighbors01 = [v1 vMidNext];
            if vMidNext ~= vMidPrev
                error('T vertex has extra neighbor!');
            end
            
            N = numNeighbors2(vMidNext);
            here = find(v1 == VV2(vMidNext,:), 1);
            VV2(vMidNext, 1:N+1) = [VV2(vMidNext, 1:here-1), vNew, ...
                VV2(vMidNext, here:N)];
            numNeighbors2(vMidNext) = N + 1;
        else
            neighbors01 = [v1 vMidNext vMidPrev];
        end
        
    end
    
    % Next handle the triangle that was bounded by (v1 v0)
    if isEdgeOnBoundary(v1, v0, VV)
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
        
        vMidNext = nextInTriangle(vNew, v0, VV2);
        vMidPrev = prevInTriangle(v1, vNew, VV2);
        
        if isOld(vMidNext) % it's a T-vertex
            neighbors10 = [v0 vMidNext];
            if vMidNext ~= vMidPrev
                error('T vertex has extra neighbor!');
            end
            
            N = numNeighbors2(vMidNext);
            here = find(v0 == VV2(vMidNext,:), 1);
            VV2(vMidNext, 1:N+1) = [VV2(vMidNext, 1:here-1), vNew, ...
                VV2(vMidNext, here:N)];
            numNeighbors2(vMidNext) = N + 1;
        else
            neighbors10 = [v0 vMidNext vMidPrev];
        end
        
    end
    
    if checkForErrors
        assert(~ismember(vNew, neighbors01));
        assert(~ismember(vNew, neighbors10));
        assert(numel(unique(neighbors01)) == numel(neighbors01));
        assert(numel(unique(neighbors10)) == numel(neighbors10));
    end
    
    % All neighbors for the vertex!!
    
    % According to the profiler, this step is verrrry slow.
    VV2(vNew, 1:(numel(neighbors01) + numel(neighbors10))) = ...
        [neighbors01 neighbors10];
    
    if checkForErrors
        if numel(unique([neighbors01 neighbors10])) <= 2
            error('poop.');
        end
    end
    
    % Now the T-matrix:
    if boundaryFlag
        T(vNew, [v0 v1]) = 0.5;
    else
        T(vNew, [v0 v1]) = 0.375;
        T(vNew, [nextInTriangle(v1, v0, VV) nextInTriangle(v0, v1, VV)]) = 0.125;
    end
    
    % Refinement matrix is always like this.
    T_refine(vNew, [v0 v1]) = 0.5;
end

%VV2_additional = sparse(newRows, newCols, newVals, ...
%    size(VV2,1), size(VV2,2));
%VV2 = VV2 + VV2_additional;

%% Tests!

getVV2 = @(v) full(VV2(v, 1:nnz(VV2(v,:))));
checkSymmetry(VV2);

Tsums = sum(T,2);

if ~all(abs(Tsums - 1.0) < 1e-6)
    error('T seems wrong.\n');
end

end


function [vertices2, VV2] = splitEdges(VV, vertices, ...
    splitAdjacentEdgeFlags, numNeighbors)
    numOriginalVertices = size(VV,1);
    vertices2 = vertices;
    VV2 = VV;

    vNew = numOriginalVertices + 1; % index of next new vertex!

    % Iterate over edges:
    %iSplit = find(splitAdjacentEdgeFlags);
    for v0 = 1:numOriginalVertices
    for jj = 1:numNeighbors(v0)
        v1 = VV(v0,jj);

        if v0 < v1 % this step prevents processing an edge twice
        if all(splitAdjacentEdgeFlags([v0 v1]))

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

            vNew = vNew + 1;
        end
        end
    end
    end
end

