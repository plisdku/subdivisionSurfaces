% Once refinement has been done, figure out which vertices belong to the
% creases and create the right transformation matrix for them.
% T_crease has size [numVertices, numOriginalVertices] and is nonzero only
% in rows corresponding to crease vertices.  Its nonzero rows should
% replace rows from subdivision of the whole mesh, for instance.
function [crease2, T_crease] = subdividedCrease(crease, ...
    creasePerturbFlags, VV, VV2)

    import VVMesh.*

    numOriginalVertices = size(VV,1);
    numVertices = size(VV2,1);
        
    assert(iscolumn(crease));

    if numel(crease) == 1 % Single fixed point
        crease2 = crease;
        T_crease = sparse(crease, crease, 1, numVertices, numOriginalVertices);

    else

        % crease(1) ~= crease(end) if it's an open-loop crease.  An
        % open-loop crease has fixed endpoints.
        if crease(1) ~= crease(end)
            creasePerturbFlags([1 end]) = false;
        end

        [crease2, T_crease] = refinedCrease(crease, creasePerturbFlags, ...
            VV, VV2);
    end

end

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
    
    numOriginalVertices = size(VV,1);
    numVertices = size(VV2,1);
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