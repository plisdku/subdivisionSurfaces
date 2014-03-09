% Loop subdivision driver


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
