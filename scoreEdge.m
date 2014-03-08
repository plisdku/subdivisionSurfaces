function key = scoreEdge(vv, ww, VV, vertices)

Nv = nnz(VV(vv,:));
Nw = nnz(VV(ww,:));

if Nv <= 3 || Nw <= 3
    key = -Inf;
    return;
end

xx = nextInTriangle(vv,ww,VV);
yy = nextInTriangle(ww,vv,VV);

key = -deltaMeanAbsoluteCurvature(vv, ww, VV, vertices);

%key = norm(vertices(ww,:) - vertices(vv,:)) / ...
%    norm(vertices(xx,:) - vertices(yy,:));

