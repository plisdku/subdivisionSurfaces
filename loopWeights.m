function omega = loopWeights(numVertices)
% omega = loopWeights(numVertices)
%
% Calculate weight factors for neighboring vertices in Loop scheme.
% These are weights for neighboring points of *existing* nodes!  For new
% nodes in mid-edge positions, weights are (3/8) for directly neighboring
% points and (1/8) for the two points distant across triangles.
%
% Weight for present node value = 1 - numVertices*omega.

omega = (5/8 - (3/8 + 0.25*cos(2*pi./numVertices)).^2)./numVertices;