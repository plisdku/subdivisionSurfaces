function [vertices, faces] = flatRegularMesh(N)
% [v f] = flatRegularMesh
% 
% Create a test mesh.

%% Create the geometry:
if nargin == 0
    N = 10;
end

vertices = zeros(N*N,3);

% Basis vectors
e1 = [1 0 0];
e2 = [0.5, 0.5*sqrt(3), 0];

% Coordinates
[ii jj] = ndgrid(1:N, 1:N);

vertices = [ii(:), jj(:)] * [e1; e2];

%% Create the topology:

idx = @(x,y) sub2ind([N N], x, y);

numFaces = (N-1)*(N-1)*2;
faces = zeros(numFaces, 3);

curFace = 1;

for y = 1:N-1
    for x = 1:N-1
        faces(curFace,:) = [idx(x,y), idx(x+1,y), idx(x,y+1)];
        faces(curFace+1,:) = [idx(x+1,y), idx(x+1,y+1), idx(x,y+1)];
        curFace = curFace + 2;
    end
end
