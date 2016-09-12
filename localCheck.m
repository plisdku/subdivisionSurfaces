function localCheck(VV, vertices)
% localCheck(VV,vertices)
%
% Check the topology of all specified vertices in VV: symmetry of adjacency
% information and whatever else I think of.  Duplicates too.

import VVMesh.*

each = @(A) reshape(A, 1, []);

for vv = each(vertices)
    
    valence = nnz(VV(vv,:));
    
    % Check valence too low
    assert(valence > 2);
    
    % Check symmetry
    for ww = each(VV(vv,1:valence))
        valence_w = nnz(VV(ww,:));
        assert(ismember(vv, VV(ww,1:valence_w)));
    end
    
    % Check for duplicates
    assert(valence == numel(unique(VV(vv,1:valence))));
    
end
