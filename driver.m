% driver!!!

%[vertices faces] = readSTL('chunk10.stl');
[vertices faces] = readSTL('chunk3.stl');
VV = fv2vv(faces(randperm(size(faces,1)),:), vertices);

figure(120); clf
patch('Faces', faces, 'Vertices', vertices, 'FaceColor', 'r')
view(3)
axis image vis3d

%% Loop subdivision!!!

[VV2 vertices2] = loopSubdivision(VV, vertices);
%%
for ii = 1:6
    %figure(ii); clf
    %plotVV(VV2, vertices2, 'b-');
    %axis image vis3d
    
    faces2 = vv2fv(VV2);
    %figure(ii+10); clf
    figure(1); clf
    pause
    fprintf('Patch time!\n');
    flatPatch('Faces', faces2, 'Vertices', vertices2, 'FaceColor', 'r');
    hold on
    plotVV(VV, vertices, 'b--')
    view(3); axis image vis3d
    camlight right
    fprintf('Pause time!\n');
    
    pause
    [VV2 vertices2] = loopSubdivision(VV2, vertices2);
end

%%
% As far as restricting certain regions to be divided or not divided:
% well... any edge that gets split can connect no matter what.  The
% splitting is a per-edge technique, that's all.  Mark edges to split that
% way!
%
% One approach then: mark a vertex as "split adjacent" or something?  Or
% split every edge next to a marked vertex?  I'm using a VV data structure
% so I gotta make this make sense for me.

%% Test mesh.

[vertices faces] = flatRegularMesh(3);
VV = fv2vv(faces, vertices);

vertices(:,3) = sin(2*(vertices(:,1) + vertices(:,2)));

%%

crease = [1:4, 5:5:20, 25:-1:22, 21:-5:1];

creaseB = [17];

%%
[VV2 vertices2 T creaseB2 crease2] = loopSubdivision(VV, vertices, creaseB, crease);
[VV2 vertices2 T2 creaseB2 crease2] = loopSubdivision(VV2, vertices2, creaseB2, crease2);
[VV2 vertices2 T3 creaseB2 crease2] = loopSubdivision(VV2, vertices2, creaseB2, crease2);
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

movObj = QTWriter('wiggle.mov');

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
    writeMovie(movObj, getframe);
    
    pause(0.01)
end
close(movObj);

%%

[VV2, ~, vertices2] = loopRefine(VV2, vertices2);
plotVV(VV2, vertices2);

%%

plotVV(VV, vertices, 'b-')
hold on
view(2); axis xy image

%%
[vertices faces] = flatRegularMesh(5);
vertices(:,3) = sin(2*(vertices(:,1) + vertices(:,2)));
VV = fv2vv(faces, vertices);

A = vv2adjacency(VV);
numVertices = size(VV,1);

refineMe = 13;
selectionFlags = zeros(numVertices,1);
selectionFlags(refineMe) = 1;

ring1 = A*selectionFlags > 0;
ring2 = A*ring1 > 0;

[VV1 v1] = truncateVV(VV, vertices, ring1);
[VV2 v2] = truncateVV(VV, vertices, ring2);

figure(1); clf
plotVV(VV, vertices, 'b-');
view(2); axis xy image

hold on
plotVV(VV2, v2, 'g-', 'LineWidth', 2)
plotVV(VV1, v1, 'r-', 'LineWidth', 3)

[VVr, ~, vr] = loopRefine(VV, vertices, ring1 | selectionFlags);
plotVV(VVr, vr, 'k-');

%% Nicer plot!

fvr = vv2fv(VVr);
figure(10); clf
patch('Faces', fvr, 'Vertices', vr, 'FaceColor', 'g', 'EdgeAlpha', 0.1)

%% Test truncation...

[vertices faces] = flatRegularMesh(3);
VV = fv2vv(faces, vertices);
keep = [6:9, 2];

[VV2 vertices2] = truncateVV(VV, vertices, keep);

%% Make a scrambled-up mesh
[vertices f] = flatRegularMesh(5);
VV = fv2vv(f,vertices);
vertices = vertices + rand(size(vertices)) - 0.5;
vertices = vertices.^2;

figure(1); clf
plotVV(VV, vertices, 'b-');
view(2);
axis image

%%