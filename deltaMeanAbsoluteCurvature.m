function deltaH = deltaMeanAbsoluteCurvature(vv, ww, VV, vertices)
% function deltaH = deltaMeanAbsoluteCurvature(vv, ww, VV, vertices)

unit = @(A) A/norm(A);

if isUnorientedEdgeOnBoundary(vv, ww, VV)
    deltaH = 0;
    return
end

xx = nextInTriangle(vv,ww,VV);
yy = nextInTriangle(ww,vv,VV);

dH_lostEdge = -calcH(vv,ww,xx,yy);
dH_gainedEdge = calcH(xx,yy,vv,ww);

dH_perimeter = calcDeltaH(vv,yy,ww,xx) + ...
    calcDeltaH(yy,ww,vv,xx) + ...
    calcDeltaH(ww,xx,vv,yy) + ...
    calcDeltaH(xx,vv,ww,yy);

deltaH = dH_lostEdge + dH_gainedEdge + dH_perimeter;

function edgeMAC = calcH(v0,v1,v2,v3)
% Calculate the MAC contribution from v0-v1, where v0-v1-v2 is a triangle
% and v3 is across from v2 (so v0-v3-v1 is a triangle too).
    
    p01 = diff(vertices([v0 v1],:));
    p02 = diff(vertices([v0 v2],:));
    p03 = diff(vertices([v0 v3],:));
    
    normal012 = unit(cross(p01, p02));
    normal031 = unit(cross(p03, p01));
    
    beta = asin(norm(cross(normal012, normal031)));
    
    edgeMAC = 0.25*abs(beta)*norm(p01);
    
    %fprintf('\tcontrib (%i %i) is %f\n', full(v0), full(v1), edgeMAC);
    
end

function deltaMAC = calcDeltaH(v0, v1, v2, v2b)
% Calculate the change of MAC along v0-v1 when its v2 is replaced with v2b.
% Will need the right-hand neighbor vertex of v0-v1!
    
    if isUnorientedEdgeOnBoundary(v0, v1, VV)
        deltaMAC = 0;
        return
    end
    
    vRHS = nextInTriangle(v1, v0, VV);
    p01 = diff(vertices([v0 v1],:));
    p02 = diff(vertices([v0 v2],:));
    p02b = diff(vertices([v0 v2b],:));
    p0R = diff(vertices([v0 vRHS],:));
    
    normal012 = unit(cross(p01, p02));
    normal012b = unit(cross(p01, p02b));
    normal0R1 = unit(cross(p0R, p01));
    
    betaOriginal = asin(norm(cross(normal0R1, normal012)));
    betaSwapped = asin(norm(cross(normal0R1, normal012b)));
    
    deltaMAC = 0.25*(abs(betaSwapped) - abs(betaOriginal))*norm(p01);
    
    %fprintf('\ttweak (%i %i) from %f to %f\n', full(v0), full(v1), ...
    %    betaOriginal, betaSwapped);
    
end

end