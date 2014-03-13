function VV2 = flipOrientation(VV, flipThese)
% VV2 = flipOrientation(VV)
%
% Reverse all vertex orderings to flip the orientation of a mesh.
%
% VV2 = flipOrientation(VV, vv)
%
% Only reverses orientation of selected vertices.
%

VV2 = VV;

if nargin == 1
    flipThese = 1:size(VV,1);
else
    flipThese = reshape(flipThese, 1, []);
end

for vv = flipThese
    neighbors = VV2(vv, 1:nnz(VV2(vv,:)));
    VV2(vv, 1:numel(neighbors)) = fliplr(neighbors);
end

