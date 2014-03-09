function [VV2 vertices2 T varargout] = loopSubdivision(VV, vertices, varargin)
% [VV2 verticesT ] = loopSubdivision(VV, vertices)
%
% Perform loop subdivision on the vertex-vertex mesh (VV, vertices)
%
% T is the transformation matrix such that vertices2 = T*vertices.
%
% [VV2 vertices2 T crease2] = loopSubdivision(VV, vertices, crease)
%
% where crease is a chain of adjacent vertices will treat the crease as a
% 1D spline in subdivision.  crease2 is the list of refined vertex
% indices making up the crease after subdivision.

import VVMesh.*

varargout = cell(numel(varargin),1);

[VV2, vertices2, T, ~, varargout{:}] = loopSubdivisionAdaptive(VV, vertices, ...
    1:size(VV,1), varargin{:});

end
