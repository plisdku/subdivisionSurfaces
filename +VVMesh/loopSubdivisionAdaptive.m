function [VV2, vertices2, T, refineVertices, varargout] = ...
    loopSubdivisionAdaptive(VV, vertices, refineVertices, varargin)
% [VV2 vertices2 T ] = loopSubdivisionAdaptive(VV, vertices, refineVertices)
%
% Perform loop subdivision on the vertex-vertex mesh (VV, vertices).
% refineVertices is a list of vertex indices bounding a set of triangles
% that must be refined.  Refinement will extend slightly outside this
% region to ensure that the inner refined region is identical to that from
% global refinement.
%
% T is the transformation matrix such that vertices2 = T*vertices.
%
% [VV2 vertices2 T refineVertices2 crease2] = 
%   loopSubdivisionAdaptive(VV, vertices, refineVertices, crease)
%
% where crease is a chain of adjacent vertices will treat the crease as a
% 1D spline in subdivision.  crease2 is the list of refined vertex
% indices making up the crease after subdivision.  refineVertices2 is the
% new list of vertices to refine if you want to continue subdivision: it
% consists of the input refineVertices as well as all new vertices bounding
% the refinement region.  It does not include the few extra subdivided
% vertices outside the refinement region.

    % TODO: rename file to loopSubdivisionLocal

    import VVMesh.*

    [VV2, ~, T, perturbFlags] = loopRefine(VV, vertices, refineVertices);
    % ok that was the easy part.

    refineVertices = union(refineVertices, ...
        neighborhood(refineVertices, VV2, 1, 2));

    %% Now do the edge splines.

    numCreases = numel(varargin);

    varargout = cell(numCreases,1);

    for cc = 1:numCreases

        crease = varargin{cc};
        [crease2, T_crease] = subdividedCrease(crease, perturbFlags(crease),...
            VV, VV2);

        T(crease2,:) = 0;
        T = T + T_crease;

        varargout{cc} = crease2;

    end

    vertices2 = T * vertices;

end





    %{
    
    assert(iscolumn(crease));
    
    if numel(crease) == 1 % Single fixed point
        crease2 = crease;
        
        T(crease,:) = 0;
        T(crease, crease) = 1.0;
    
    else
        creasePerturbFlags = perturbFlags(crease);
        
        % crease(1) ~= crease(end) if it's an open-loop crease.  An
        % open-loop crease has fixed endpoints.
        if crease(1) ~= crease(end)
            creasePerturbFlags([1 end]) = false;
        end
        
        [crease2, T_crease] = refinedCrease(crease, creasePerturbFlags, ...
            VV, VV2); 
        
        T(crease2,:) = 0;
        T = T + T_crease;
    end
    
    varargout{cc} = crease2;
    %}
