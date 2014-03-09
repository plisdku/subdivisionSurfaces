function plotVV(VV, vertices, varargin)

numVertices = size(vertices, 1);
numEdges = nnz(VV)/2;


pStarts = zeros(numEdges, size(vertices,2));
pEnds = zeros(numEdges, size(vertices,2));

iEdge = 1;

for ii = 1:numVertices
    %fprintf('v %i\n', ii);
    jj = VV(ii, VV(ii,:) > ii); % Only take each edge once.
    
    nPlotEdges = numel(jj);
    
    % Plot edge (ii,jj) with jj > ii.
    pStarts(iEdge:iEdge+nPlotEdges-1,:) = vertices(repmat(ii,1,nPlotEdges),:);
    pEnds(iEdge:iEdge+nPlotEdges-1,:) = vertices(jj,:);
    
    iEdge = iEdge + nPlotEdges;
end

%% Now plot!

heldAlready = ishold();

numDims = size(vertices,2);

if numDims == 2
    
    plot([pStarts(:,1) pEnds(:,1)]', ...
        [pStarts(:,2) pEnds(:,2)]', ...
        varargin{:});
    hold on
    %plot3(vertices(:,1), vertices(:,2), 'o');
    
elseif numDims == 3
    
    plot3([pStarts(:,1) pEnds(:,1)]', ...
        [pStarts(:,2) pEnds(:,2)]', ...
        [pStarts(:,3) pEnds(:,3)]', ...
        varargin{:});
    hold on
    %plot3(vertices(:,1), vertices(:,2), vertices(:,3), 'o');
    
else
    error('Num dims = %i, weirdo.', numDims);
end

if false == heldAlready
    hold off
end