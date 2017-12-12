function plotFV(faces, vertices, varargin)

vv = VVMesh.fv2vv(faces, vertices);
VVMesh.plotVV(vv, vertices, varargin{:});