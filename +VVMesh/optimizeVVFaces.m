function VV2 = optimizeVVFaces(VV, vertices, fixedEdges)
% VV2 = optimizeVVFaces(VV, vertices)
%
% Iteratively flip diagonals of split quads to make triangles more uniform.

import VVMesh.*

checkNoDuplicates(VV);
checkSymmetry(VV);
        
numVertices = size(VV,1);

valence = numVVNeighbors(VV);
numEdges = sum(valence)/2;

makePlots = false;
checkForErrors = false;

%% Create matrix of fixed edges

numFixedEdges = size(fixedEdges,1);
fixed = sparse(fixedEdges(:,1), fixedEdges(:,2), true(numFixedEdges,1),...
    numVertices, numVertices);
fixed = fixed | transpose(fixed);

%% Prioritize edges for swapping.  (Build priority queue.)
% Put ALL non-boundary edges into a priority queue.
% 
% Edges connected to vertices of valence 3 (or less) are deemed unswappable
% but may become swappable through the course of optimization (if their
% valence grows larger than 3).  So put those edges into the priority queue
% with key -Inf.
%
% Other edges should be given a key based on how swappable they are, e.g.
% based on curvature reduction or triangle uniformity.

edges = zeros(numEdges, 2);
scores = zeros(numEdges, 1);

currEdge = 1;

for vv = 1:numVertices
    
    for ii = 1:valence(vv)
    ww = VV(vv,ii);
    if vv < ww && ~isUnorientedEdgeOnBoundary(vv,ww,VV) && ~fixed(vv,ww)
        
        score = scoreEdge(vv, ww, VV, vertices);
        
        edges(currEdge,:) = [vv ww];
        scores(currEdge,:) = score;
        
        currEdge = currEdge + 1;
    end
    end
    
end

nEdge = currEdge-1;

%% Make the priority queue: a heap

heap = Heap(scores(1:nEdge), edges(1:nEdge,:));

%% While the top score > 1, keep swapping edges.

VV2 = VV;

if heap.heapSize == 0
    return;
end

done = false;

while ~done
    
    [score, edge] = heap.top();
    
    if score <= 1
        done = true;
    else
        % The two existing triangles are (vv, ww, xx) and (ww, vv, yy)
        
        vv = edge(1);
        ww = edge(2);
        Nv = nnz(VV2(vv,:)); % valence of v
        Nw = nnz(VV2(ww,:)); % valence of w
        
        xx = nextInTriangle(vv,ww,VV2);
        yy = nextInTriangle(ww,vv,VV2);
        Nx = nnz(VV2(xx,:));
        Ny = nnz(VV2(yy,:));
        
        %figure(2); clf
        %plotSomeVV(VV2, vertices, [vv ww xx yy], 2, 'b-');
        %hold on
        %plotSomeVV(VV2, vertices, [vv ww xx yy], 'r-', 'LineWidth', 2);
        %plotSomeVV(VV2, vertices, [vv ww], 'g--', 'LineWidth', 3);
        
        if checkForErrors
            assert(~ismember(xx, VV(yy, 1:Ny)));
            assert(~ismember(yy, VV(xx, 1:Nx)));
            assert(Nv > 3);
            assert(Nw > 3);
            localCheck(VV2, [vv ww xx yy]);
        end
        
        %fprintf('Processing %f, presently %f, going to %f\n', score, ...
        %    norm(diff(vertices([vv ww],:))),  ...
        %    norm(diff(vertices([xx yy],:))));
        
        % Plot first:
        
        if makePlots
            f2 = vv2fv(VV2);

            figure(1); clf
            patch('Faces', f2, 'Vertices', vertices, 'FaceColor', 'g', ...
                'EdgeAlpha', 0.1);
            camlight right
            lighting phong
            %view(3);
            view(90,-40);
            axis image
            [truncVV truncVertices] = truncateVV(VV2, vertices, [xx yy vv ww]);
            hold on
            plotVV(truncVV, truncVertices, 'r-');
            %plotVV(trVV, trv, 'y--', 'LineWidth', 4)
            pause(0.4);
        end
        
        % vv removes ww, ww removes vv
        
        VV2(vv, 1:Nv) = [setdiff(VV2(vv,1:Nv),ww, 'stable') 0];
        %checkNoDuplicates(VV2);
        VV2(ww, 1:Nw) = [setdiff(VV2(ww,1:Nw),vv, 'stable') 0];
        %checkNoDuplicates(VV2);
        
        % x adds y in between v and w
        
        v_in_x = find(VV2(xx,1:Nx) == vv, 1);
        VV2(xx, 1:Nx+1) = [VV2(xx, 1:v_in_x) yy VV2(xx, v_in_x+1:Nx)];
        
        % y adds x in between w and v
        
        w_in_y = find(VV2(yy,1:Ny) == ww, 1);
        VV2(yy, 1:Ny+1) = [VV2(yy, 1:w_in_y) xx VV2(yy, w_in_y+1:Ny)];
        
        
        % Check that new edge has vertices of valence > 3
        if checkForErrors
            localCheck(VV2, [vv ww xx yy]);
        
            assert(nnz(VV2(xx,:)) > 3);
            assert(nnz(VV2(yy,:)) > 3);
        end
        
        % Rescore it: score -> 1/score
        
        heap.pop();
        heap.push(1/score, sort([xx yy]));
        
        if checkForErrors
            assert(ww == nextInTriangle(xx,yy,VV2));
            assert(vv == nextInTriangle(yy,xx,VV2));
        end
        
        % Function to score an edge:
        keyFunction = @(edge) scoreEdge(VV2, vertices, edge(1), edge(2));
        
        % Rescore the 1-neighborhood of the swapped edge: all edges from
        % (xx or yy) to (vv or ww).
        
        % Rescore the 2-neighborhood of the swapped edge (2-nhood of new
        % edge).  This consists of edges from (vv or ww) to a point
        % neighboring (xx or yy), and edges from (xx or yy) to a point
        % neighboring (vv or ww).  Get these edges:
        
        % Get all edges from p1 to a neighbor of p2.
        % This function is gross but effective.
        makeNeighborhoodEdges = @(p1, p2) cell2mat(...
            arrayfun(@(p2) [p1 p2],...
                full(intersect(nonzeros(VV2(p1,:)), nonzeros(VV2(p1,:)))), ...
                'UniformOutput', false));
        
        N2EdgesUnsorted = [makeNeighborhoodEdges(vv,xx); ...
            makeNeighborhoodEdges(vv,yy); ...
            makeNeighborhoodEdges(ww,xx); ...
            makeNeighborhoodEdges(ww,yy)];
        N2Edges = unique(sort(N2EdgesUnsorted,2), 'rows');
        
        for ee = 1:size(N2Edges,1)
            %fprintf('Rey-keying [%i %i] to %f\n', ...
            %    N2Edges(ee,1), N2Edges(ee,2), keyFunction(N2Edges(ee,:)));
            
            heap.reKey(N2Edges(ee,:), ...
                scoreEdge(N2Edges(ee,1), N2Edges(ee,2), VV2, vertices));
        end
        
    end
    
end
