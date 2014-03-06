% driver!!!

%[vertices faces] = readSTL('chunk10.stl');
[vertices faces] = readSTL('chunk3.stl');

figure(120); clf
patch('Faces', faces, 'Vertices', vertices, 'FaceColor', 'r')
view(3)
axis image vis3d

%%

VV = fv2vv(faces(randperm(size(faces,1)),:), vertices);
%VV = fv2vv(faces, vertices);

%% Test selective splitting

splitFlags = vertices(:,3) < -50;
%splitFlags = vertices(:,3) < 250;
A = vv2adjacency(VV);
adjFlags = (A^1)*splitFlags > 0;

%%

[VV2, ~, vertices2] = loopRefine(VV, vertices, adjFlags);
f2 = vv2fv(VV2);


%%
figure(1); clf
plotVV(VV2, vertices2, 'b-')

%%
figure(120); clf
subplot(121)
plotVV(VV, vertices, 'b-')
hold on
plot3(vertices(adjFlags,1), vertices(adjFlags,2), vertices(adjFlags,3),...
    'ro', 'LineWidth', 3)
view(3); axis image vis3d

subplot(122)
plotVV(VV2, vertices2, 'b-')
hold on
plotVV(VV, vertices, '-', 'Color', [0.8 0.8 1])
view(3); axis image vis3d
%%

%%

figure(9); clf
plotVV(VV, vertices, 'b-', 'LineWidth', 2)
[VV2, ~, vertices2] = loopRefine(VV, vertices);
[VV2, ~, vertices2] = loopRefine(VV2, vertices2);
[VV2, ~, vertices2] = loopRefine(VV2, vertices2);

hold on
plotVV(VV2, vertices2, 'r--')

%%

faces2 = vv2fv(VV);

figure(20); clf
patch('Faces', faces2, 'Vertices', vertices, 'FaceColor', 'g')
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
    figure(ii+10); clf
    flatPatch('Faces', faces2, 'Vertices', vertices2, 'FaceColor', 'r');
    hold on
    plotVV(VV, vertices, 'b--')
    view(3); axis image vis3d
    camlight right
    
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

[vertices faces] = flatRegularMesh(5);
VV = fv2vv(faces, vertices);

vertices(:,3) = sin(2*(vertices(:,1) + vertices(:,2)));

%%

[VV2 vertices2] = loopSubdivision(VV, vertices);

figure(1); clf
plotVV(VV, vertices, 'b-');
hold on
plotVV(VV2, vertices2, 'r-');
view(2)
axis image vis3d

%%
[VV2, ~, vertices2] = loopRefine(VV2, vertices2);
plotVV(VV2, vertices2);

