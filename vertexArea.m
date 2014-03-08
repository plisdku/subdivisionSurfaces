function area = vertexArea(vv, VV, vertices)
% function area = vertexArea(vv, VV, vertices)

area = 0;

neighbors = VV(vv, 1:nnz(VV(vv,:)));
numNeighbors = numel(neighbors);

wrap = @(n) 1 + mod(n-1, numNeighbors);

p0 = vertices(vv,:);

for ii = 1:numNeighbors
    jj = wrap(ii+1);
    
    if ~isEdgeOnBoundary(vv, neighbors(ii), VV)

        p1 = vertices(neighbors(ii),:);
        p2 = vertices(neighbors(jj),:);

        area = area + norm(cross(p1-p0, p2-p0)); % quad area! not tri area.
        
        %fprintf('area = %f\n', area);
    end
    
end

% The area is the sum of areas of parallelograms.  Divide by two to get
% areas of triangles.  Divide by three to distribute areas evenly across
% all vertices of each triangle.  Net result: must divide by 6.

area = area/6;

