% Driver for mesh optimization

%% Make a distorted mesh and optimize it
[vertices f] = flatRegularMesh(10);
VV = fv2vv(f,vertices);
vertices(:,3) = 4*(sin(2*(vertices(:,1) + vertices(:,2))));
vertices = vertices + 0.25*(rand(size(vertices)) - 0.5);
vertices = vertices.^2;

[VV vertices] = loopSubdivision(VV, vertices);
VV2 = optimizeVVFaces(VV, vertices);

%%
figure(1); clf
plotVV(VV, vertices, 'b-');
view(2);
%axis normal
axis image

figure(2); clf
plotVV(VV2, vertices, 'b-');
view(2);
%axis normal
axis image

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

%% Find new triangles only.

% Step 1: get a unique representation of faces.

uf1 = faces;
uf2 = f2;
numFaces = size(faces,1);

for ff = 1:numFaces
    uf1(ff,:) = sort(uf1(ff,:));
    uf2(ff,:) = sort(uf2(ff,:));
end

% Step 2: find the changed faces.
[uniqueIn1, if1] = setdiff(uf1, uf2, 'rows');
[uniqueIn2, if2] = setdiff(uf2, uf1, 'rows');

figure(3); clf
%plotVV(VV, vertices, 'Color', [0.8 0.8 1]);
hold on
flatPatch('Vertices', vertices, 'Faces', faces(uniqueIn1,:), ...
    'FaceColor', 'r', 'FaceAlpha', 0.25, 'EdgeAlpha', 0.1);
flatPatch('Vertices', vertices, 'Faces', f2(uniqueIn2,:), ...
    'FaceColor', 'g', 'FaceAlpha', 0.25, 'EdgeAlpha', 0.1);
axis xy image vis3d
view(3)

%%

nTheta = 100;

thetas = linspace(0, 360, nTheta+1);

while 1
for ii = 1:nTheta
    [az el] = view;
    view([az+thetas(2)-thetas(1), el]);
    pause(0.01)
    
end
end

%%
figure(3); clf
plotVV(VV, vertices, 'Color', [0.8 0.8 1]);
hold on
plotVV(VV2, vertices, 'Color', [1 0.1 0.1])

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