% Driver for local refinement.

import VVMesh.*

[vertices faces] = flatRegularMesh(5);
vertices(:,3) = sin(2*(vertices(:,1) + vertices(:,2)));
VV = fv2vv(faces, vertices);

refineVertices = [1 2 6];

% Refinement means:
% - split all edges in the 1-neighborhood of the refinement vertices
% - the new refinement vertices are in the 1-neighborhood of AT LEAST TWO
% of the old refinement vertices.  Use the adjacency matrix to figure this
% out, since it can count the connections between vertices.

%%

VV2 = VV;
vertices2 = vertices;

A = vv2adjacency(VV2);
    
Ttot = [];
numRefinements = 10;
for refinements = 1:numRefinements
    
    % Refine all edges between vertices in the 1-neighborhood of the
    % refinement set.
    
    [VV3, vertices3, T, ~, T_refine] = loopRefine(VV2, vertices2, refineVertices);
    
    % Test ...
    v3_matrixVersion = T_refine*vertices2;
    assert(max(abs(vertices3(:) - v3_matrixVersion(:))) < 1e-12);
    
    % Perform Loop interpolation
    vertices3 = T*vertices2;
    refineVertices = union(neighborhood(refineVertices, VV3, 1, 2), ...
        refineVertices);
    
    % The new vertices to refine are the original vertices to refine and
    % all the new (split) vertices on the edges between them.
    
    
    VV2 = VV3;
    vertices2 = vertices3;
    
    figure(1); clf
    plotVV(VV, vertices, 'b-'); %, 'Color', [0.8 0.8 1]);
    hold on
    plotVV(VV2, vertices2, 'b-');
    plotSomeVV(VV2, vertices2, refineVertices, 'r-');
    axis xy image vis3d
    view(2)
    pause
end