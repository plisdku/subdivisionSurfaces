function [VV2 v2 T indexKey varargout] = simplify(VV, vertices, costFunction, varargin)
% [VV2 v2 T key] = simplify(VV, vertices, costFunction)
% [VV2 v2 T key] = simplify(VV, vertices, minimumLength)
% [VV2 v2 T key crease2] = simplify(VV, vertices, costFunction, crease)
%
% Collapse edges to simplify the mesh.
%
% T is the matrix such that v2 = T*vertices.
%
% key is a lookup table to tell you what happened to all the vertices.
%
% costFunction should have the signature @(vv,ww,VV,vertices).
% The edge with the largest cost function will be collapsed first.
% Once all edges have cost functions <= 0, simplification stops.
% 
%

import VVMesh.*
makePlots = false;
checkForErrors = false;

if makePlots
    figNum = gcf;
end

if nargin < 3 || isempty(costFunction)
    minLength = 5;
    costFunction = @(vv, ww, VV, verts) ...
        minLength^2 - sum( (verts(vv,:)-verts(ww,:)).^2 );
    assert(~isempty(costFunction(1,1,VV,vertices)));
elseif isnumeric(costFunction)
    minLength = costFunction;
    costFunction = @(vv, ww, VV, verts) ...
        minLength^2 - sum( (verts(vv,:)-verts(ww,:)).^2 );
    assert(~isempty(costFunction(1,1,VV,vertices)));
end

if nargin < 4
    fixedEdges = [];
else
    fixedEdges = creaseEdges(varargin);
end

if checkForErrors
    checkNoDuplicates(VV);
    checkSymmetry(VV);
    assert(isOrientedConsistently(VV));
end
        
numVertices = size(VV,1);

valence = numVVNeighbors(VV);
numEdges = sum(valence)/2;


%% Create matrix of fixed edges

numFixedEdges = size(fixedEdges,1);

if numFixedEdges > 0
    fixed = sparse(fixedEdges(:,1), fixedEdges(:,2), true(numFixedEdges,1),...
        numVertices, numVertices);
    fixed = fixed | transpose(fixed);
else
    fixed = sparse(numVertices, numVertices);
end

fixedVertexFlags = false(numVertices,1);
fixedVertexFlags(fixedEdges(:)) = true;

%% Prioritize edges for swapping.  (Build priority queue.)
% Put ALL non-boundary edges into a priority queue.
% Edges are scored according to the cost function.

edges = zeros(numEdges, 2);
scores = zeros(numEdges, 1);

currEdge = 1;

for vv = 1:numVertices
    
    for ii = 1:valence(vv)
    ww = VV(vv,ii);
    if vv < ww && ~fixed(vv,ww)
        
        score = costFunction(vv, ww, VV, vertices);
        
        edges(currEdge,:) = [vv ww];
        scores(currEdge,:) = score;
        
        currEdge = currEdge + 1;
    end
    end
    
end

nEdge = currEdge-1;

%fprintf('Num scores > 0: %i\n', nnz(scores > 0));

%% Make the priority queue: a heap

heap = Heap(scores(1:nEdge), edges(1:nEdge,:));

% This is how I can keep track of what happens to all those edges!
ordering = heap.swapBuffer(1:heap.numSwaps);

% Look up where something is in the heap.  Very low-level and hacky.  :-p
% heapIdx(ii,jj) is the heap array index of the edge from ii to jj.
heapIdx = sparse(heap.values(:,1), heap.values(:,2), 1:heap.heapSize, ...
    numVertices, numVertices);

deletedEdges = spalloc(numVertices,numVertices,heap.heapSize);

%% Discarded vertices end up here:

deletedVerts = zeros(numVertices,1);
numDeleted = 0;

%% Matrix to create new vertices

T = speye(numVertices);

%% While the top score > 0, keep swapping edges.
p3 = @(v, varargin) plot3(v(:,1), v(:,2), v(:,3), varargin{:});

VV2 = VV;
v2 = vertices;

if heap.heapSize == 0
    return;
end

done = false;


while ~done
    
    [score, edge] = heap.top();
    
    if deletedEdges(edge(1),edge(2))
        %fprintf('Completing lazy deletion from heap.\n');
        
        heap.pop();
        updateHeapIndex();
        
        if checkForErrors, checkHeapIndex(); end
    elseif nearValenceThree(edge)
        
        %fprintf('Keeping edge neighboring a vertex of valence 3\n');
        
        heap.pop();
        updateHeapIndex();
        
        if checkForErrors, checkHeapIndex(); end
    elseif all(fixedVertexFlags(edge))
        %fprintf('Edge has become fixed because both ends are fixed.\n');
        heap.pop();
        updateHeapIndex();
        if checkForErrors, checkHeapIndex(); end
    elseif score <= 0
        done = true;
    else
        score2 = costFunction(edge(1), edge(2), VV2, v2);
        assert(score == score2);
        % Collapse works by removing vertex ww and keeping vertex vv.
        % If ww is fixed, though, so it can't be moved, then I should
        % delete vv.  Or, swap vv and ww!  That's what I do.  :-)
        if fixedVertexFlags(edge(2))
            assert(~fixedVertexFlags(edge(1)));
            vv = edge(2);
            ww = edge(1);
        else
            vv = edge(1);
            ww = edge(2);
        end
        
        %if any(edge == 107)
        %    fprintf('Got you.\n');
        %end
        
        Nv = nnz(VV2(vv,:)); % valence of v
        Nw = nnz(VV2(ww,:)); % valence of w
        
        %fprintf('%i to %i\n', vv, ww);
        
        vNeigh = VV2(vv,1:Nv);
        wNeigh = VV2(ww,1:Nw);
        
        if makePlots
            % PLOT 1:
            figure(figNum);
            subplot(221); cla
            plotSomeVV(VV2, v2, neighborhood([vv;ww], VV2), 'bo-');
            hold on; p3(v2([vv;ww],:), 'ro', 'LineWidth', 4);
            %plotVVLabels(VV2, v2, neighborhood([vv;ww],VV2), ...
            %    'FontSize', 18);
            axis vis3d; view(3)
            ax = axis;
            ax1 = gca;
            pause(0.001)
        end
        
        % 1. Add ww to the orphans list.
        numDeleted = numDeleted + 1;
        deletedVerts(numDeleted) = ww;
        
        % 2. For all neighbors of w: if v is not present already, replace w
        % with v.  (Since v and w have become subsumed by v, I only need to
        % keep one instance of v.)
        %
        % Inside the heap, the edge ending on w must be re-valued to end on
        % v instead.
        
        for neighborVert = setdiff(wNeigh, vv)
            oldEdge = sort([neighborVert ww]);
            if ~any(VV2(neighborVert,:) == vv)
                %fprintf('Replacing %i in %i\n', full(ww), full(neighborVert));
                VV2(neighborVert, VV2(neighborVert,:) == ww) = vv;
                
                newEdge = sort([neighborVert vv]);
                
                oldIdx = heapIdx(oldEdge(1), oldEdge(2));
                heapIdx(newEdge(1), newEdge(2)) = oldIdx;
                heap.values(oldIdx,:) = newEdge;
            end
            %heapIdx(oldEdge(1), oldEdge(2)) = DELETED;
            deletedEdges(oldEdge(1), oldEdge(2)) = 1;
        end
        
        % 3. For all neighbors of v: if w is present, delete it.
        % This amounts to deleting a bunch of edges from the heap, in fact.
        % Somehow this has to be dealt with!
        
        for neighborVert = setdiff(vNeigh, ww)
            iw = find(VV2(neighborVert,:) == ww);
            if iw
                Nn = nnz(VV2(neighborVert,:));
                VV2(neighborVert, 1:Nn) = ...
                    [VV2(neighborVert, [1:iw-1, iw+1:Nn]), 0];
                
                oldEdge = sort([neighborVert ww]);
                %heapIdx(oldEdge(1), oldEdge(2)) = DELETED;
                deletedEdges(oldEdge(1), oldEdge(2)) = 1;
            end
        end
        
        % 3. Insert ww's neighbors into vv, where ww was.
        % We do this in VV2.  The heap is just waiting to be re-keyed.
        
        v_in_w = find(VV2(ww,1:Nw) == vv, 1);
        w_in_v = find(VV2(vv,1:Nv) == ww, 1);
        
        newVNeigh = [ VV2(vv, 1:w_in_v-1), ...
            VV2(ww, [v_in_w+1:Nw, 1:v_in_w-1]), ...
            VV2(vv, w_in_v+1:Nv) ];
        newVNeigh = unique(newVNeigh, 'stable');
        VV2(vv, 1:numel(newVNeigh)) = newVNeigh;
        VV2(vv, numel(newVNeigh)+1:end) = 0;
        assert(nnz(VV2(vv,:)) == numel(newVNeigh));
        
        % 4.  Clear the ww row from VV2.  I need to keep all deleted rows
        % until the end to avoid horrible array re-sizing and re-indexing
        % runtimes.  Clearing rows will keep the adjacency matrix
        % symmetrical and leave the deleted points just dangling in space,
        % but present on the deleted vertex list so I can handle this issue
        % later.
        VV2(ww,:) = 0;
        
        % 5.  Move v2(vv,:) to a better position!
        pv0 = v2(vv,:); % save for plotting later.
        if ~fixedVertexFlags(vv)
            v2(vv,:) = 0.5*(v2(vv,:) + v2(ww,:));
            
            T(vv,:) = 0.5*T(vv,:) + 0.5*T(ww,:);
        end
        
        if checkForErrors
            assert(~any(VV2(:) == ww));
            localCheck(VV2, vv);
            localCheck(VV2, neighborhood(vv, VV2));
            checkSymmetry(VV2);
            checkNoDuplicates(VV2);
            assert(isOrientedConsistently(VV2));
        end
        
        if makePlots
            figure(figNum);
            subplot(222); cla
            plotSomeVV(VV2, v2, [vv; neighborhood(vv, VV2)], 'bo-');
            hold on; p3(v2(vv,:), 'ro', 'LineWidth', 4);
            %plotVVLabels(VV2, v2, [vv; neighborhood(vv,VV2)], ...
            %    'FontSize', 18);
            axis vis3d; view(3)
            axis(ax);
            ax2 = gca;
            
            subplot(2,2,[3 4]); cla
            plotVV(VV2, v2, 'bo-');
            hold on; p3(v2(vv,:), 'ro', 'LineWidth', 4);
            p3([pv0; v2(ww,:)], 'rx', 'LineWidth', 4);
            %plotVVLabels(VV2, v2, [vv; neighborhood(vv,VV2)], ...
            %    'FontSize', 18);
            axis image vis3d; view(3)
            %axis(ax)
            
            %linkprop([ax1 ax2], {'CameraPosition', 'CameraTarget', ...
            %    'CameraUpVector', 'CameraViewAngle'});
            
            pause
            %pause(0.001)
        end
        
        if checkForErrors, checkHeapIndex(); end
        
        heap.pop();
        deletedEdges(vv,ww) = 1;
        updateHeapIndex();
        
        if checkForErrors, checkHeapIndex(); end
        
        % Now re-key the heap!
        % All edges attached to vv must change.
        
        vNeigh = nonzeros(VV2(vv,:))';
        for nv = vNeigh
            
            newEdge = sort([vv nv]);
            
            newKey = costFunction(newEdge(1), newEdge(2), VV2, v2);
            heap.reKeyIndex(heapIdx(newEdge(1), newEdge(2)), newKey);
            updateHeapIndex();
            
            %fprintf('Performed %i swaps in re-key.\n', heap.numSwaps);
        end
        
        if checkForErrors, checkScores(); end
        
        if makePlots
            pause(0.1)
        end
        
        
    end
    
    if heap.heapSize == 0
        %fprintf('Heap is empty: we are done.\n');
        done = true;
    end
    
end

% Eliminate deleted verts!

deletedVerts = deletedVerts(1:numDeleted);
numNewVertices = numVertices - numDeleted;
keptVerts = setdiff(1:numVertices, deletedVerts);
indexKey = zeros(numVertices,1);
indexKey(keptVerts) = 1:numNewVertices;
indexKey(deletedVerts) = numNewVertices+1:numVertices;

% Push all the deleted vertices to the end, everywhere.
% First, translate all indices in VV2 and translate row numbers.
[rr cc ss] = find(VV2);
VV2 = sparse(indexKey(rr), cc, indexKey(ss), numVertices, numVertices);

if checkForErrors
    checkNoDuplicates(VV2);
    checkSymmetry(VV2);
    assert(isOrientedConsistently(VV2));
end

% The transformation matrix gets a row-reordering but the columns must stay
% the same to line up with the input vector.
[rr cc ss] = find(T);
T = sparse(indexKey(rr), cc, ss, numVertices, numVertices);

% Next, re-order the vertices:
[~,iReorder] = sort(indexKey);
v2 = v2(iReorder,:);

% Finally truncate everything!
VV2 = VV2(1:numNewVertices,:);
T = T(1:numNewVertices,:);
v2 = v2(1:numNewVertices,:);

if checkForErrors
    checkNoDuplicates(VV2);
    checkSymmetry(VV2);
    assert(isOrientedConsistently(VV2));
end


% Translate crease indices when they've been disrupted by the vertex
% consolidation.
varargout = cell(size(varargin));

for cc = 1:numel(varargin)
    newCrease = indexKey(varargin{cc});
    varargout{cc} = newCrease(newCrease <= numNewVertices);
end

return;





% nested:
function updateHeapIndex()
    
    % This can actually be quite fast!  I don't need to do any swaps, I can
    % just read everything straight out of the heap.  Thank goodness.
    
%    heapIdx = sparse(heap.values(1:heap.heapSize,1),...
%        heap.values(1:heap.heapSize,2), ...
%        1:heap.heapSize, ...
%        numVertices, numVertices);
    
    % List of all heap entries I need to go adjust.  I don't seem to need
    % to know the order that things were swapped in because all the
    % information I need is in heap.values() at the locations in
    % touchList.
    %touchList = unique(heap.swapBuffer(:));
    touchList = heap.swapBuffer(1:heap.numSwaps,:); % can be redundant, it won't hurt!!
    
    p1 = heap.values(touchList(:),1);
    p2 = heap.values(touchList(:),2);
    
    if checkForErrors, assert(all(p1 < p2)); end % edge ordering rule!
    
    % The effect here is, for EACH p1 and p2,
    %  heapIdx(p1(n),p2(n)) = touchList(n).
    % The only exception should be deleted edges!
    nii = sub2ind(size(heapIdx), p1, p2);
    %nii = p1 + (p2-1)*size(heapIdx,1);
    
    %nonDeleted = heapIdx(nii) ~= DELETED;
    %heapIdx(nii(nonDeleted)) = touchList(nonDeleted);
    
    heapIdx(nii) = touchList;
    
    
end

function checkHeapIndex()
    % Probably the best way to do this is to run through the heap and make
    % sure that heapIdx contains the right information.  heapIdx also has
    % out-of-date records for edges that had vertices collapsed and such,
    % and those records don't need to be checked.
    
    for nii = 1:heap.heapSize
        edge = heap.values(nii,:);
        
        % The heap can contain edges that have been deleted from the mesh.
        % They are lazily removed from the heap as they bubble to the top.
        % Deleted edges are marked as DELETED in the heapIdx matrix, and at
        % least one of their two vertices will then be in the deleted
        % vertex list.
        %if heapIdx(edge(1),edge(2)) == DELETED
        if deletedEdges(edge(1), edge(2))
            assert(any(ismember(edge, deletedVerts)))
        else
            assert(heapIdx(edge(1),edge(2)) == nii);
        end
    end
end

% Determine whether the vertices in edge are two out of three neighbors of
% any adjacent vertex.  If so then collapsing the edge will leave that
% vertex with valence two, which is forbidden.
function yesNo = nearValenceThree(edge)
    
    sharedNeighbors = intersect(nonzeros(VV2(edge(1),:)), ...
        nonzeros(VV2(edge(2),:)));
    
    valences = arrayfun(@(ww) nnz(VV2(ww,:)), sharedNeighbors);
    
    yesNo = any(valences <= 3);
    
end

function checkScores
    for hh = 1:heap.heapSize
        myEdge = heap.values(hh,:);
        
        if ~deletedEdges(myEdge(1), myEdge(2))
            myScore = heap.keys(hh,:);
            myScore2 = costFunction(myEdge(1), myEdge(2), VV2, v2);

            assert(myScore == myScore2);
        end
    end
    
    fprintf('All scores ok.\n');
end

end



function edges = creaseEdges(creases)
    
    edges = [];
    
    col = @(A) reshape(A, [], 1);
    
    for cc = 1:length(creases)
        edges = vertcat(edges,...
            [col(creases{cc}(1:end-1)) col(creases{cc}(2:end))]);
    end
end