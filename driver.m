% driver!!!

import VVMesh.*

myPatch = @(f,v,color) flatPatch('Vertices', v, 'Faces', f, ...
    'FaceColor', color, 'EdgeAlpha', 0.1);
p3 = @(v, varargin) plot3(v(:,1), v(:,2), v(:,3), varargin{:});

[vertices faces] = readSTL('chunk3.stl');

VV = fv2vv(faces, vertices);
numVertices = size(VV,1);

figure(1); clf
plotSomeVV(VV, vertices, find(vertices(:,3) < 0), 'b-');
axis image
view(10, 60)

%%

topVerts = find(vertices(:,3) == 0);
botVerts = find(vertices(:,3) == -150);

topCrease = topVerts(boundaryVertices(truncateVV(VV, topVerts)));
botCrease = botVerts(boundaryVertices(truncateVV(VV, botVerts)));

figure(2); clf
plotSomeVV(VV, vertices, topVerts, 'Color', [0.8 0.8 1]);
hold on
%plot3(vertices(topVerts,1), vertices(topVerts,2), vertices(topVerts,3), ...
%    'LineWidth', 2)
plot3(vertices(topCrease,1), vertices(topCrease,2), vertices(topCrease,3),...
    'LineWidth', 2)

%% Now that I've got the top and bottom creases, fix and subdivide.

[VV2 v2 T top2 bot2] = loopSubdivision(VV, vertices, ...
    topCrease([1:end,1]), botCrease([1:end,1]));

%%

f2 = vv2fv(VV2);
figure(1); clf
myPatch(f2, v2, 'r')
hold on
p3(v2(top2,:), 'b-', 'LineWidth', 3)
p3(v2(bot2,:), 'b-', 'LineWidth', 3)

%plotVV(VV2, v2, 'b-')

%% Do a subdivision that doesn't touch the top!

topNeighborhood = neighborhood(topVerts, VV);
refineVerts = setdiff(1:numVertices, topNeighborhood);

[VV2 v2 T rv2 top2 bot2] = loopSubdivisionAdaptive(VV, vertices, refineVerts,...
    topCrease([1:end,1]), botCrease([1:end,1]));
[VV2 v2 T rv2 top2 bot2] = loopSubdivisionAdaptive(VV2, v2, rv2,...
    top2, bot2);


