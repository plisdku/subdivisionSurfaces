function alpha = edgeEdgeAngles(uu, VV, vertices)

valence = nnz(VV(uu,:));
alpha = zeros(valence,1);
wrap = @(n) 1 + mod(n-1, valence);

for ee = 1:valence
    
    if ~isEdgeOnBoundary(uu, VV(uu,ee), VV)
        v0 = vertices(VV(uu,ee),:) - vertices(uu,:);
        v1 = vertices(VV(uu,wrap(ee+1)),:) - vertices(uu,:);
        alpha(ee) = acos(dot(v0,v1)/norm(v0)/norm(v1));
    else
        alpha(ee) = NaN;
    end
    
end

alpha = alpha(~isnan(alpha));