% Test curvature functions!

import VVMesh.*

testCloseTol = @(a, b, tol) assert(norm(a-b) < tol);
testClose = @(a, b) testCloseTol(a, b, 1e-6);

%% Flat wagon-wheel mesh

[vertices faces] = wagonWheel(4);
VV = fv2vv(faces, vertices);

%% Vertex area!

% Central vertex area: 2/3
testClose(vertexArea(5, VV, vertices), 2/3);

% Edge vertices: area 1/3
for ii = 1:4
    testClose(vertexArea(ii, VV, vertices), 1/3);
end

fprintf('Area tests PASSED!\n');

%% Normal vectors!

testClose(normalVector(1,2,VV,vertices), [0 0 1]);

try
    normalVector(2,1,VV,vertices);
    error('Should have failed!');
catch err; end

testClose(normalVector(5,1,VV,vertices), [0 0 1]);

fprintf('Easy normal vector tests PASSED!\n');

%% Angles!

alphas = edgeEdgeAngles(5,VV,vertices);
testClose(alphas, [pi/2, pi/2, pi/2, pi/2]');

alphas = edgeEdgeAngles(1,VV,vertices);
testClose(alphas, [pi/4, pi/4]');

testClose(dihedralAngle(5,1,VV,vertices), 0);

fprintf('Easy angle tests PASSED!\n');

%% Curvatures!
% The Gaussian curvature is not defined (?) for edge vertices.
testClose(vertexGaussianCurvature(5, VV, vertices), 0);

testClose(vertexMeanCurvature(5, VV, vertices), 0);

fprintf('Easy curvature tests PASSED!\n');

%% Now test some things on a more angled mesh.
% Use the corner of a cube:

unit = @(v) v/norm(v);

[~, faces] = wagonWheel(3);
vertices = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
VV = fv2vv(faces, vertices);

testClose(normalVector(4,1,VV,vertices), [0 0 1]);
testClose(normalVector(4,2,VV,vertices), [1 0 0]);
testClose(normalVector(4,3,VV,vertices), [0 1 0]);

fprintf('Hard normal vector tests PASSED!\n');

%%

alphas = edgeEdgeAngles(4, VV, vertices);
testClose(alphas, [pi/2 pi/2 pi/2]');

testClose(dihedralAngle(4,1,VV,vertices), pi/2);

fprintf('Hard angle tests PASSED!\n');

%%

testClose(vertexGaussianCurvature(4, VV, vertices), (2*pi - 3*pi/2));

testClose(vertexMeanCurvature(4, VV, vertices), 0.25*3*pi/2);

fprintf('Hard curvature tests PASSED!\n');

%% Test total mean absolute curvature perturbation!

[vertices faces] = flatRegularMesh(5);
numVertices = size(vertices,1);
vertices(:,3) = rand(numVertices,1);

VV = fv2vv(faces,vertices);

Htot = totalMeanAbsoluteCurvature(VV, vertices);

%% Now flip the edges and see what happens.

for vv = 1:numVertices
    Nv = nnz(VV(vv,:));
    for ee = 1:Nv
    	ww = VV(vv,ee);
        Nw = nnz(VV(ww,:)); % valence of w
        
        % MAKE A COPY of the topology.
        VV2 = VV;
        
        if isUnorientedEdgeOnBoundary(vv,ww,VV) || Nv <= 3 || Nw <= 3 || vv > ww
            continue;
        end
        
        xx = nextInTriangle(vv,ww,VV2);
        yy = nextInTriangle(ww,vv,VV2);
        Nx = nnz(VV2(xx,:));
        Ny = nnz(VV2(yy,:));
        
        %figure(1); clf
        %plotVV(VV2, vertices, 'b-')
        %view(2); axis image;
        %hold on
        %plotSomeVV(VV2, vertices, [xx yy ww vv], 'r-', 'LineWidth', 2);
        
        assert(~ismember(xx, VV(yy, 1:Ny)));
        assert(~ismember(yy, VV(xx, 1:Nx)));
        
        % vv removes ww, ww removes vv
        VV2(vv, 1:Nv) = [setdiff(VV2(vv,1:Nv),ww, 'stable') 0];
        %checkNoDuplicates(VV2);
        VV2(ww, 1:Nw) = [setdiff(VV2(ww,1:Nw),vv, 'stable') 0];
        %checkNoDuplicates(VV2);
        
        
        
        % x adds y in between v and w
        v_in_x = find(VV2(xx,1:Nx) == vv, 1);
        VV2(xx, 1:Nx+1) = [VV2(xx, 1:v_in_x) yy VV2(xx, v_in_x+1:Nx)];
        
        % y adds x in between w and v
        w_in_y = find(VV2(yy,1:Ny) == ww, 1);
        VV2(yy, 1:Ny+1) = [VV2(yy, 1:w_in_y) xx VV2(yy, w_in_y+1:Ny)];
        
        % Checks
        
        %assert(ww == nextInTriangle(xx,yy,VV2));
        %assert(vv == nextInTriangle(yy,xx,VV2));
        %localCheck(VV2, [vv ww xx yy]);
        
        Htot2 = totalMeanAbsoluteCurvature(VV2, vertices);
        
        deltaH_expected = deltaMeanAbsoluteCurvature(vv, ww, VV, vertices);
        deltaH_obtained = Htot2 - Htot;
        
        testClose(deltaH_expected, deltaH_obtained);
        
        %fprintf('Delta = %f (expected %f)\n', deltaH_obtained, ...
        %    deltaH_expected);
        
    end
end


fprintf('Delta total mean curvature test PASSED!\n');





