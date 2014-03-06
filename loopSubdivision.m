function [VV2 vertices2 T varargout] = loopSubdivision(VV, vertices, varargin)
% [VV2 verticesT ] = loopSubdivision(VV, vertices)
%
% Perform loop subdivision on the vertex-vertex mesh (VV, vertices)
%
% T is the transformation matrix such that vertices2 = T*vertices.
%
% [VV2 vertices2 T crease2] = loopSubdivision(VV, vertices, crease)
%
% where crease is a chain of adjacent vertices will treat the crease as a
% 1D spline in subdivision.  crease2 is the list of refined vertex
% indices making up the crease after subdivision.

numOriginalVertices = size(vertices,1);

%if nargin < 3
%    splitAdjacentEdgeFlags = ones(numOriginalVertices, 1);
%end

numNeighbors = numVVNeighbors(VV);

[VV2, T, ~, inheritance] = loopRefine(VV, vertices);
% ok that was the easy part.

numVertices = size(VV2,1);

%% Now do the edge splines.

numCreases = numel(varargin);

for cc = 1:numCreases
    
    crease = varargin{cc};
    crease2 = refinedCrease(crease, VV, VV2);
    
    % Wipe out T for the crease
    T(crease2,:) = 0;
    
    % Set T for the new vertices
    Tnew = sparse(crease2([2:2:end, 2:2:end]), crease([1:end-1, 2:end]),...
        0.5, numVertices, numOriginalVertices);
    
    % Set T for the old vertices
    
    if numel(crease) == 1 % hard point
        Told = sparse(crease, crease, 1, numVertices, numOriginalVertices);
    elseif crease(1) == crease(end) % boundary loop
        Told = sparse(crease([1:end-1,1:end-1,1:end-1]), ...
            crease([end-1, 1:end-2, 1:end-1, 2:end-1, 1]), ...
            [0.125*ones(numel(crease)-1,1), ...
             0.75*ones(numel(crease)-1,1), ...
             0.125*ones(numel(crease)-1,1)], ...
             numVertices, numOriginalVertices);
    else
        Told_ends = sparse(crease([1 end]), crease([1 end]), [1 1], ...
            numVertices, numOriginalVertices);
        
        Told_mid = sparse(crease([2:end-1, 2:end-1, 2:end-1]), ...
            crease([1:end-2, 2:end-1, 3:end]), ...
            [0.125*ones(numel(crease)-2,1), ...
             0.75*ones(numel(crease)-2,1), ...
             0.125*ones(numel(crease)-2,1)], ...
             numVertices, numOriginalVertices);
        
        Told = Told_ends + Told_mid;
    end
    
    T = T + Tnew + Told;
    
    varargout{cc} = crease2;
end

%%

vertices2 = T * vertices;

end


function refined = refinedCrease(crease, VV, VV2)
% Get list of all vertices in the crease after refinement.
%
% crease is the list of vertices in the crease before refinement.
%
% VV and VV2 are the topologies pre- and post- refinement.

findIndex = @(iFindMe, indices) find(indices == iFindMe, 1);

refined = zeros(2*numel(crease) - 1, 1);
refined(1:2:end) = crease(:);

for ii = 2:2:numel(refined)
    
    v0 = refined(ii-1);
    v1 = refined(ii+1);
    
    vMid = VV2(v0, findIndex(v1, VV(v0,:)));
    assert(isscalar(vMid));
    
    refined(ii) = vMid;
    
end


end
