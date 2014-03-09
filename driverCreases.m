% Driver for subdivision with creases.

%% Test mesh.

[vertices faces] = flatRegularMesh(5);
VV = fv2vv(faces, vertices);

vertices(:,3) = sin(2*(vertices(:,1) + vertices(:,2)));

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
