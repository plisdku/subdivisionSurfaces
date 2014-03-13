% test orientation!

import VVMesh.*

[vertices faces] = flatRegularMesh(4);
VV = fv2vv(faces);

% Vertices are in a staggered array, connected to nearest neighbors:
%
%       13  14  15  16
%     9   10  11  12
%   5   6   7   8
% 1   2   3   4
%

%% Nonadjacent vertices cause errors (unless valence 2)
% I'm not sure whether I would like valence 2 or the error condition to
% take precedence.  Technically, valence 2 wins right?

assert(~threw(@() areSameOrientation(1,16,VV))); % valence 2, no errors
assert(threw(@() areSameOrientation(2, 7, VV))); % non-adjacent, > 2
assert(threw(@() areSameOrientation(7, 2, VV))); % non-adjacent, > 2
assert(~threw(@() areSameOrientation(2,3, VV))); % adjacent: no errrors

fprintf('Errors on nonadjacent vertices PASSED!\n');

%% Vertices with valence < 3 have both orientations

assert(areSameOrientation(1, 2, VV)); % true because v1 has valence 2
assert(areSameOrientation(1, 5, VV)); % true because v1 has valence 2
assert(areSameOrientation(5, 1, VV)); % symmetric

VV2 = flipOrientation(VV, 1);

assert(areSameOrientation(1, 2, VV2)); % true because v1 has valence 2
assert(areSameOrientation(1, 5, VV2)); % true because v1 has valence 2
assert(areSameOrientation(5, 1, VV2)); % symmetric

fprintf('Valence < 3 PASSED!\n');

%% Adjacent vertices of valence > 2

assert(areSameOrientation(2, 3, VV)); % two edge vertices
assert(areSameOrientation(3, 2, VV));

assert(areSameOrientation(2, 6, VV)); % one edge vertex, one inner vertex
assert(areSameOrientation(6, 2, VV));

assert(areSameOrientation(10, 11, VV)); % two inner vertices
assert(areSameOrientation(11, 10, VV));

VV2 = flipOrientation(VV, [3 6 11]);

assert(~areSameOrientation(2, 3, VV2)); % two edge vertices
assert(~areSameOrientation(3, 2, VV2));

assert(~areSameOrientation(2, 6, VV2)); % one edge one inner
assert(~areSameOrientation(6, 2, VV2)); % one edge one inner 

assert(areSameOrientation(3, 6, VV2)); % two flipped vertices
assert(areSameOrientation(6, 3, VV2));

assert(~areSameOrientation(10, 11, VV2)); % two inner vertices
assert(~areSameOrientation(11, 10, VV2));

fprintf('Valence >= 3 PASSED!\n');

%% Check that mesh orientation flipper touches and flips the right vertices

flipThese = [9 10 8]';
VV2 = flipOrientation(VV, flipThese);

[VVFlipped, flipped, touched] = orient(VV2, 6);

assert(isequal(sort(find(flipped)), sort(flipThese)));

valence = numVVNeighbors(VV);
assert(isequal(sort(find(valence > 2)), sort(find(touched))));
assert(isequal(sort(find(valence < 3)), sort(find(~touched))));

fprintf('Flipping and touching test PASSED!\n');

% And of course check that the whole mesh is appropriately oriented.

for vv = 1:size(VV,1)
    for ww = nonzeros(VV(vv,:))'
        assert(areSameOrientation(vv, ww, VVFlipped));
    end
end

fprintf('Re-oriented mesh orientation PASSED!\n');

assert(isOrientedConsistently(VV));
assert(~isOrientedConsistently(VV2));
assert(isOrientedConsistently(VVFlipped));

fprintf('isOrientedConsistently tests PASSED!\n');
