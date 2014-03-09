% driver!!!

%N = 15;
%[vertices faces] = flatRegularMesh(N);
%vertices(:,3) = sin(0.07*vertices(:,1).^2);

[vertices faces] = readSTL('chunk3.stl');

VV = fv2vv(faces, vertices);
numVertices = size(VV,1);

%vertices(:,3) = sin(2*(0.5*vertices(:,1).^2 + vertices(:,2)));


figure(1); clf
plotSomeVV(VV, vertices, find(vertices(:,3) < 0), 'b-');
axis image
view(10, 60)

%%
[trVV trVertices] = truncateVV(VV, vertices, find(vertices(:,3) < 0));
trFaces = vv2fv(trVV);

figure(1); clf
patch('Faces', trFaces, 'Vertices', trVertices, 'FaceColor', 'g', ...
    'EdgeAlpha', 0.1);
axis image vis3d
view(10, 60)
camlight left

%%

selection = num2cell(find(vertices(:,3) < -60));
[VV2 vertices2] = loopSubdivision(VV, vertices, selection{:});

%%


[trVV trv] = truncateVV(VV, vertices, find(vertices(:,3) < 0));

Hv = zeros(size(trVV,1), 1);
for vv = 1:size(trVV,1)
    %Hv(vv) = vertexMeanCurvature(vv, VV, vertices);
    Hv(vv) = vertexMeanCurvature(vv, trVV, trv);
end

%%

trFaces = vv2fv(trVV);

figure(2); clf
%patch('Faces', faces, 'Vertices', vertices, 'FaceVertexCData', Hv, ...
patch('Faces', trFaces, 'Vertices', trv, 'FaceVertexCData', Hv, ...
    'FaceColor', 'interp');
colorbar
axis image vis3d
view(10, 60)

%%

vertices2 = vertices;
VV2 = VV;

refineVertices = find(vertices2(:,1) > -20 & ...
    vertices2(:,1) < 60 & ...
    vertices2(:,2) > -20 & ...
    vertices2(:,2) < 20 & ...
    vertices2(:,3) > -40 & ...
    vertices2(:,3) < 0);
rv2 = refineVertices;

for itr = 1:3
    
    [VV2, vertices2, ~, rv2] = loopSubdivisionAdaptive(VV2, vertices2, ...
        rv2);
    
    [trVV trVertices] = truncateVV(VV2, vertices2,  ...
        find(vertices2(:,3) < 0));
    trFaces = vv2fv(trVV);

    figure(1); clf
    patch('Faces', trFaces, 'Vertices', trVertices, 'FaceColor', 'g', ...
        'EdgeAlpha', 0.1);
    axis image vis3d
    view(10, 60)
    camlight left
    lighting phong
    zoom(3)
    hold on
    plot3(vertices2(rv2,1), vertices2(rv2,2), vertices2(rv2,3), 'ro');
    
    pause

end
%%
%refineVertices = find(Hv > 30);

figure(3); clf
plotVV(VV2, v2, 'b-');

%%


%%

Hv2 = zeros(size(VV2,1), 1);
for vv = 1:size(VV2,1)
    Hv2(vv) = vertexMeanCurvature(vv, VV2, v2);
end

%%

f2 = vv2fv(VV2);

figure(102); clf
patch('Faces', f2, 'Vertices', v2, 'FaceVertexCData', Hv2, ...
    'FaceColor', 'interp');
colorbar
axis image vis3d
view(10, 60)

%%

crease = [1:4, 5:5:20, 25:-1:22, 21:-5:1]';

rv0 = neighborhood([8 12 13], VV, 2);

%%
[VV2 vertices2 T rv crease2] = loopSubdivisionAdaptive(VV, vertices, rv0, crease);
%%
[VV2 vertices2 T2 rv crease2] = loopSubdivisionAdaptive(VV2, vertices2, rv, crease2);
%%
[VV2 vertices2 T3 rv crease2] = loopSubdivisionAdaptive(VV2, vertices2, rv, crease2);
%%
f2 = vv2fv(VV2);

figure(1); clf
whitebg k
plotVV(VV, vertices, 'b-');
hold on
patch('Faces', f2, 'Vertices', vertices2, 'FaceColor', 'r', ...
    'FaceAlpha', 1.0, 'EdgeAlpha', 0.1);
%plotVV(VV2, vertices2, 'r-');
plot3(vertices(crease,1), vertices(crease,2), vertices(crease,3), ...
    'b--', 'LineWidth', 3)
plot3(vertices2(crease2,1), vertices2(crease2,2), vertices2(crease2,3),...
    'g-', 'LineWidth', 4)
plot3(vertices2(creaseB2,1), vertices2(creaseB2,2), vertices2(creaseB2,3),...
    'g-', 'LineWidth', 4)
axis image vis3d
axis off
set(gcf, 'Color', 'k')
view(3);
camlight right; lighting phong
ax = axis;

%%

Ttot = T3*T2*T;
verts = vertices;

%movObj = QTWriter('wiggle.mov');

for tt = 1:100
    
    %randVert = randi(25,1);
    %verts(randVert,:) = verts(randVert,:) + rand - 0.5;
    
    verts = verts + 0.1*rand(size(verts)) - 0.05;
    verts2 = Ttot*verts;
    
    figure(10); clf
    plotVV(VV, verts, 'b-');
    hold on
    patch('Faces', f2, 'Vertices', Ttot*verts, 'FaceColor', 'r', ...
        'EdgeAlpha', 0.1);
    plot3(verts(crease,1), verts(crease,2), verts(crease,3), ...
        'b--', 'LineWidth', 3)
    plot3(verts2(crease2,1), verts2(crease2,2), verts2(crease2,3),...
        'g-', 'LineWidth', 4)
    plot3(verts2(creaseB2,1), verts2(creaseB2,2), verts2(creaseB2,3),...
        'g-', 'LineWidth', 4)
    axis off
    axis(ax)
    view(az, el)
    set(gcf, 'Color', 'k')
    camlight right; lighting phong
    %writeMovie(movObj, getframe);
    
    pause(0.01)
end
%close(movObj);

