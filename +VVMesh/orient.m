function [VV2 flipped touched] = orient(VV, v0)
% VV2 = orient(VV, v0)
%
% Propagate the orientation of vertex v0 across the mesh VV in an attempt
% to make it consistently oriented.

    import VVMesh.*

    assert(isscalar(v0));
    numVertices = size(VV,1);
    valence = numVVNeighbors(VV);
    
    if valence(v0) < 3
        error('Cannot orient to a vertex of valence < 3.');
    end
    
    touched = false(numVertices,1);
    touched(v0) = true; % this means, don't visit it.
    flipped = false(numVertices,1);
    
    VV2 = VV;
    
    stack = zeros(20, 2); % [oriented, neighborIndex]
    stackSize = 0;
    
    push(v0, 1);
    
    while ~emptyStack()
        
        [vCurr ii] = pop();
        ww = VV2(vCurr,ii);
        
        % If vCurr may pass its orientation to ww:
        if ww ~= 0
            push(vCurr, ii+1);
            if ~touched(ww) && valence(ww) > 2
                if ~areSameOrientation(vCurr, ww, VV2)
                    VV2(ww, 1:valence(ww)) = fliplr(VV2(ww, 1:valence(ww)));
                    flipped(ww) = true;
                end
                
                touched(ww) = true;
                push(ww, 1);
            end
        end
        
    end
    

    function tf = emptyStack()
        tf = ~logical(stackSize);
    end
    
    function push(vCurr, idx)
        stackSize = stackSize + 1;
        
        if length(stack) < stackSize % dynamically resize stack RARELY
            actualSize = size(stack, 1);
            newSize = max(2*actualSize, numVertices);
            stack = [stack; zeros(newSize - actualSize, 2)];
        end
        
        stack(stackSize,:) = [vCurr idx];
    end

    function [vCurr ii] = pop()
        vCurr = stack(stackSize,1);
        ii = stack(stackSize, 2);
        stackSize = stackSize - 1;
    end

end
