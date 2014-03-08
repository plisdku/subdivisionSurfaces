function [vertices faces] = wagonWheel(nOuter)
% function [vertices faces] = wagonWheel(nOuter)

if nargin == 0
    nOuter = 3;
end

numVertices = nOuter + 1;

vertices = zeros(numVertices, 3);

thetas = linspace(0, 2*pi, numVertices);

for nn = 1:nOuter
    vertices(nn,1:2) = [cos(thetas(nn)) sin(thetas(nn))];
end

faces = zeros(nOuter,3);

for tt = 1:nOuter
    wrap = @(n) 1 + mod(n-1, nOuter);
    faces(tt,:) = [numVertices, wrap(tt), wrap(tt+1)];
end
