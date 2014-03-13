function [VV2 vertices2 T refineVertices varargout] = loopSubdivisionAdaptive(...
    VV, vertices, refineVertices, varargin)
% [VV2 verticesT ] = loopSubdivisionAdaptive(VV, vertices, refineVertices)
%
% Perform loop subdivision on the vertex-vertex mesh (VV, vertices)
%
% T is the transformation matrix such that vertices2 = T*vertices.
%
% [VV2 vertices2 T refineVertices2 crease2] = 
%   loopSubdivisionAdaptive(VV, vertices, refineVertices, crease)
%
% where crease is a chain of adjacent vertices will treat the crease as a
% 1D spline in subdivision.  crease2 is the list of refined vertex
% indices making up the crease after subdivision.

import VVMesh.*

numOriginalVertices = size(vertices,1);

[VV2, T, vertices2, perturbFlags] = loopRefine(VV, vertices, refineVertices);
% ok that was the easy part.

refineVertices = union(refineVertices, ...
    neighborhood(refineVertices, VV2, 1, 2));

numVertices = size(VV2,1);

%% Now do the edge splines.

numCreases = numel(varargin);

varargout = cell(numCreases,1);

for cc = 1:numCreases
    
    crease = varargin{cc};
    
    assert(iscolumn(crease));
    
    if numel(crease) == 1 % Single fixed point
        crease2 = crease;
        
        T(crease,:) = 0;
        T(crease, crease) = 1.0;
    
    elseif crease(1) == crease(end) % Boundary loop
        creasePerturbFlags = perturbFlags(crease);
        % the boundary loop might be fixed externally at the "endpoints",
        % so don't set the ends to be movable!!
        %creasePerturbFlags([1 end]) = true;
        
        [crease2 T_crease] = refinedCrease(crease, creasePerturbFlags, ...
            VV, VV2);
        
        T(crease2,:) = 0;
        T = T + T_crease;
    else
        creasePerturbFlags = perturbFlags(crease);
        
        % only fixed ends make sense.  i think.
        creasePerturbFlags([1 end]) = false;
        
        [crease2 T_crease] = refinedCrease(crease, creasePerturbFlags, ...
            VV, VV2);
        
        T(crease2,:) = 0;
        T = T + T_crease;
    end
    
    varargout{cc} = crease2;
end

%vertices2 = T * vertices;



function [refined, T] = refinedCrease(crease, perturbFlags, VV, VV2)
    % For closed boundary loops:
    % call with crease(1) == crease(end), and perturbFlags == 1 at ends.
    %
    % For fixed-end boundaries:
    % call with crease(1) ~= crease(end), and perturbFlags == 0 at ends.
    
    if crease(1) == crease(end) % boundary loop
        creaseLength = numel(crease) - 1;
        assert(perturbFlags(1) == perturbFlags(end));
    else
        creaseLength = numel(crease);
        assert(~any(perturbFlags([1 end])));
    end
    
    refined = zeros(2*numel(crease), 1); % allocate max possible length.
    
    %% Create T matrix for the original vertices.
    % Fixed vertices stay fixed; vertices that may move take a weighted
    % average of neighbors and self.
    
    wrap = @(n) 1 + mod(n-1, creaseLength);
    
    % Indices into crease of perturbed and fixed vertices.
    % Using 1:creaseLength will discard the redundant last element in the
    % case that we're handling a boundary loop.
    iPerturb = find(perturbFlags(1:creaseLength));
    iFixed = find(~perturbFlags(1:creaseLength));
    
    numPerturb = numel(iPerturb);
    vPerturb = crease(iPerturb);
    vPrev = crease(wrap(iPerturb - 1));
    vNext = crease(wrap(iPerturb + 1));
    
    numFixed = numel(iFixed);
    vFixed = crease(iFixed);
    ii = [vPerturb; vPerturb; vPerturb; vFixed];
    jj = [vPrev; vPerturb; vNext; vFixed];
    ss = [ones(numPerturb,1)/8; ...
        ones(numPerturb,1)*3/4; ...
        ones(numPerturb,1)/8; ...
        ones(numFixed,1)];
    
    T = sparse(ii, jj, ss, numVertices, numOriginalVertices);
    
    
    assert(max(sum(T,2)) <= 1.00001);
    assert(nnz(T(numOriginalVertices+1:end,:)) == 0);
    
    %% Insert midpoints into the refined crease array.
    % Also put averaging weights into T.  This is a little slow but who
    % cares.
    
    iDst = 1; % index into refined crease: destination index
    
    % Check it out... numel(crease)-1 is the last vertex to put a midpoint
    % after, regardless of whether it's a loop or not.
    for iSrc = 1:numel(crease)-1

        v0 = crease(iSrc);
        v1 = crease(iSrc + 1);
        
        % If v0 and v1 are not neighbors then there is a midpoint!
        
        if ismember(v0, VV2(v1,:))
            assert(ismember(v1, VV2(v0,:)));
            
            refined(iDst) = v0;
            iDst = iDst + 1;
        else
            % We'll attempt to find a midpoint.  If there isn't one, then
            % we take the crease straight through without a midpoint.
            %
            % The midpoint vertex should be the only new vertex common to
            % v0 and v1.
        
            vMid = intersect(setdiff(VV2(v0,:), VV(v0,:)), ...
                setdiff(VV2(v1,:), VV(v1,:)));
            
            assert(numel(vMid) == 1);
            assert(vMid > numOriginalVertices);
            
            %fprintf('vMid = %i, adding %i and %i\n', full(vMid), full(v0), full(v1));
            
            refined([iDst, iDst+1]) = [v0 vMid];
            iDst = iDst + 2;
            
            % Crease matrix for NEW vertices: an average.
            T(vMid, v0) = 1/2;
            assert(max(sum(T,2)) <= 1.00001);
            T(vMid, v1) = 1/2;
            assert(max(sum(T,2)) <= 1.00001);
        end

    end

    %% Terminate the new crease.
    % If it's a boundary loop, then crease(end) == crease(1) is
    % appropriate.  If it's not a closed path then I need the end too.
    
    refined(iDst) = crease(end);
    refined = refined(1:iDst); % chop off extra length if needed

    assert(max(sum(T,2)) <= 1.00001);
end

end