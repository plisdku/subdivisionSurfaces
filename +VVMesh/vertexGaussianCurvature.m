function Kv = vertexGaussianCurvature(uu, VV, vertices)
% function K = vertexGaussianCurvature(uu, VV, vertices)
%
% Calculate the integral Gaussian curvature in the area attributed to
% vertex uu.

import VVMesh.*

Kv = 2*pi - sum(edgeEdgeAngles(uu, VV, vertices));

%Kv = K / vertexArea(uu, VV, vertices);