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

refineVertices = [1 2 6];

vertices2 = vertices;
VV2 = VV;

% Refinement means:
% - split all edges in the 1-neighborhood of the refinement vertices
% - the new refinement vertices are in the 1-neighborhood of AT LEAST TWO
% of the old refinement vertices.  Use the adjacency matrix to figure this
% out, since it can count the connections between vertices.

%%

A = vv2adjacency(VV2);
    
Ttot = [];
numRefinements = 10;
for refinements = 1:numRefinements
    
    % Refine all edges between vertices in the 1-neighborhood of the
    % refinement set.
    
    [VV2, T] = loopRefine(VV2, vertices2, refineVertices);
    vertices2 = T*vertices2;
    refineVertices = union(neighborhood(refineVertices, VV2, 1, 2), ...
        refineVertices);
    
    % The new vertices to refine are the original vertices to refine and
    % all the new (split) vertices on the edges between them.
    
    
    figure(1); clf
    plotVV(VV, vertices, 'b-'); %, 'Color', [0.8 0.8 1]);
    hold on
    plotVV(VV2, vertices2, 'b-');
    plotSomeVV(VV2, vertices2, refineVertices, 'r-');
    axis xy image vis3d
    view(2)
    pause
end
%%

[VV1 v1] = truncateVV(VV, vertices, ring1);
[VV2 v2] = truncateVV(VV, vertices, ring2);

figure(1); clf
plotVV(VV, vertices, 'b-');
view(2); axis xy image

hold on
plotVV(VV2, v2, 'g-', 'LineWidth', 2)
plotVV(VV1, v1, 'r-', 'LineWidth', 3)

[VVr, T, vr] = loopRefine(VV, vertices, ring1 | refineVertexFlags);
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

%% Make a distorted mesh and optimize it
[vertices f] = flatRegularMesh(10);
VV = fv2vv(f,vertices);
vertices(:,3) = 4*(sin(2*(vertices(:,1) + vertices(:,2))));
vertices = vertices + 0.25*(rand(size(vertices)) - 0.5);
vertices = vertices.^2;

[VV vertices] = loopSubdivision(VV, vertices);

figure(1); clf
plotVV(VV, vertices, 'b-');
view(2);
axis normal

VV2 = optimizeVVFaces(VV, vertices);

figure(2); clf
plotVV(VV2, vertices, 'b-');
view(2);
axis normal

%% Optimize a more difficult mesh

[vertices faces] = readSTL('chunk3.stl');
VV = fv2vv(faces, vertices);

%%

[VV vertices] = loopSubdivision(VV, vertices);
faces = vv2fv(VV);

%%
VV2 = optimizeVVFaces(VV, vertices);

%%

H1 = totalMeanAbsoluteCurvature(VV, vertices);
H2 = totalMeanAbsoluteCurvature(VV2, vertices);
fprintf('Total mean absolute curvature went from %f to %f\n', H1, H2);

%%
f2 = vv2fv(VV2);

figure(1); clf
flatPatch('Faces', faces, 'Vertices', vertices, 'FaceColor', 'r');
axis image vis3d;
view(3)
camlight right

figure(2); clf
flatPatch('Faces', f2, 'Vertices', vertices, 'FaceColor', 'r');
axis image vis3d;
view(3)
camlight right

%% Plot the per-vertex curvature

numVertices = size(VV,1);
Hv = zeros(numVertices,1);
Kv = Hv;

for vv = 1:numVertices
    Hv(vv) = vertexMeanCurvature(vv,VV,vertices);
    Kv(vv) = vertexGaussianCurvature(vv,VV,vertices);
end

%%

vColor = linspace(0, 1, numVertices)' * [1 1 1];

figure(3); clf
patch('Faces', faces, 'Vertices', vertices, 'FaceVertexCData', Hv, ...
    'FaceColor', 'interp', 'EdgeAlpha', 0);
axis image vis3d;
view(3);
colormap hot
colorbar
%camlight right


%% spokes

[vertices faces] = wagonWheel(4);
VV = fv2vv(faces, vertices);

VV2 = optimizeVVFaces(VV, vertices);





%%