function bVerts = boundaryVertices(VV, vSelected)
% ww = boundaryVertices(VV)
%
% Finds a boundary vertex in VV and then chains vertices counterclockwise
% until it returns to the starting vertex.
%
% ww = boundaryVertices(VV, vSelected)
%
% Finds a boundary chain in the subset of VV given by indices vSelected.
%

    import VVMesh.*

    if nargin == 2
        bVerts = vSelected(boundaryOfEntireMesh(...
            truncateVV(VV, vSelected)));
    else
        bVerts = boundaryOfEntireMesh(VV);
    end

end




function bVerts = boundaryOfEntireMesh(trVV)

    import VVMesh.*
    
    numVertices = size(trVV, 1);

    found = false;

    for vv = 1:numVertices
        for ww = nonzeros(trVV(vv,:))'
            if isEdgeOnBoundary(vv,ww,trVV)
                found = true;
                break
            end
        end  
        if found
            break
        end
    end

    if ~found
        error('Could not find any boundary edges!');
    end

    bVerts = zeros(numVertices,1);
    bVerts(1) = vv;

    prevVert = vv;
    currVert = ww;
    ii = 2;

    while currVert ~= vv
        nextVert = nextInTriangle(prevVert, currVert, trVV);

        bVerts(ii) = currVert;

        ii = ii + 1;
        prevVert = currVert;
        currVert = nextVert;
    end

    bVerts = bVerts(1:ii-1);
    
end
