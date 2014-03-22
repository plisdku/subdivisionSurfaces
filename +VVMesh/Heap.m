classdef Heap < handle
    
    properties
        keys = [];
        heapSize = 0;
        values = [];
        isRow = true;
        
        % Things that change frequently...
        
        numSwaps = 0;
        swapBuffer = [];
    end
    
    methods 
        
        function obj = Heap(varargin)
            
            if nargin == 0
                % do nothing
            else
                [obj.keys ii] = sort(varargin{1}, 'descend');
                obj.heapSize = numel(obj.keys);
                
                obj.clearSwapBuffer();
                obj.numSwaps = obj.heapSize;
                obj.swapBuffer = [ (1:obj.numSwaps)' reshape(ii, [], 1) ];
                
                if isrow(obj.keys)
                    obj.isRow = true;
                else
                    obj.isRow = false;
                end
                
                if numel(varargin) > 1
                    obj.values = varargin{2};
                    
                    if obj.isRow
                        obj.values = obj.values(:,ii);
                    else
                        obj.values = obj.values(ii,:);
                    end
                end
            end
            
        end
        
        function clearSwapBuffer(obj)
            obj.numSwaps = 0;
            
            if numel(obj.swapBuffer) < numel(obj.keys)
                obj.swapBuffer(numel(obj.keys),2) = 0;
            end
        end
        
        % T is the heap reordering matrix.
        function aOut = swap(obj, a, b)
            
            obj.keys([a b]) = obj.keys([b a]);
            
            if ~isempty(obj.values)
                if obj.isRow
                    obj.values(:,[a b]) = obj.values(:, [b a]);
                else
                    obj.values([a b],:) = obj.values([b a],:);
                end
            end
            
            obj.numSwaps = obj.numSwaps + 1;
            obj.swapBuffer(obj.numSwaps,:) = [a b];
            
            aOut = b;
            
        end
        
        function v = value(obj, idx)
            
            if obj.isRow
                v = obj.values(:,idx);
            else
                v = obj.values(idx,:);
            end
            
        end
        
        function push(obj, key, value)
            
            obj.heapSize = obj.heapSize + 1;
            iLast = obj.heapSize;
            
            % Efficient heap re-sizing.
            if obj.heapSize > numel(obj.keys) && obj.heapSize > 1
                
                newAllocSize = 2*numel(obj.keys);
                obj.keys(newAllocSize) = 0; % double size!
                
                if nargin == 3
                    if obj.isRow
                        obj.values(end,newAllocSize) = 0;
                    else
                        obj.values(newAllocSize,end) = 0;
                    end
                end     
            end
            
            obj.keys(iLast) = key;
            
            if nargin == 3
                if obj.heapSize == 1
                    % If we're inserting the first element, it's too early
                    % to know if the arrays grow in the row direction or
                    % the column direction.  If value is a column vector,
                    % then it should accumulate in the row direction.
                    obj.isRow = ~isrow(value);
                end
                
                if obj.isRow
                    obj.values(:,iLast) = value;
                else
                    obj.values(iLast,:) = value;
                end
            end
            
            obj.clearSwapBuffer();
            obj.bubbleUp(iLast);
            
        end
        
        function [key val] = top(obj)
            
            if obj.heapSize > 0
                key = obj.keys(1);
                
                if ~isempty(obj.values)
                    if obj.isRow
                        val = obj.values(:,1);
                    else
                        val = obj.values(1,:);
                    end
                end
            else
                error('Heap is empty.');
            end
            
        end
        
        function [key val] = pop(obj)
            if obj.heapSize < 1
                error('Heap is empty.');
            end
            
            if nargout > 1
                [key val] = obj.top();
            else
                key = obj.top();
            end
            
            obj.clearSwapBuffer();
            obj.swap(1, obj.heapSize);
            obj.heapSize = obj.heapSize - 1;
            
            obj.bubbleDown(1);
            
        end
        
        % Pull the key-value pair at heap index iCurrent as far up the heap
        % as it can go.  (Re-heapify when a key is larger than its parent
        % key.)
        function bubbleUp(obj, iCurrent)
            parent = @(n) floor(n/2);
            while iCurrent > 1 && ...
                obj.keys(iCurrent) > obj.keys(parent(iCurrent))

                iCurrent = swap(obj, iCurrent, parent(iCurrent));
            end
        end
        
        % Push the key-value pair at heap index iCurrent as far down the
        % heap as it can go.  (Re-heapify when a key is smaller than one of
        % its child keys.)
        function bubbleDown(obj, iCurrent)
            
            leftChild = @(n) 2*n;
            rightChild = @(n) 2*n + 1;
            
            done = false;
            while ~done
                %lc = leftChild(iCurrent);
                %rc = rightChild(iCurrent);
                
                lc = 2*iCurrent;
                rc = lc + 1;
                
                currKey = obj.keys(iCurrent);
                leftKey = -Inf;
                rightKey = -Inf;
                
                if lc <= obj.heapSize
                    leftKey = obj.keys(lc);
                end
                
                if rc <= obj.heapSize
                    rightKey = obj.keys(rc);
                end
                
                if rightKey > leftKey
                    if rightKey > currKey
                        iCurrent = swap(obj, iCurrent, rc);
                    else
                        done = true;
                    end
                else
                    if leftKey > currKey
                        iCurrent = swap(obj, iCurrent, lc);
                    else
                        done = true;
                    end
                end
            end
        end
        
        % Search the entire value array for a row that == valueAsRow
        function idx = findRow(obj, valueAsRow)
            
            idx = find(obj.values(:,1) == valueAsRow(1));
            
            if numel(idx) > 1
                dims = numel(valueAsRow);
                for dd = 1:dims
                    idx = idx(obj.values(idx,dd) == valueAsRow(dd));
                    
                    if numel(idx) == 1
                        return;
                    end
                end
            end
            
            if numel(idx) > 1
                warning('Found multiple.');
                idx = idx(1);
            end
            
        end
        
        % Search the entire value array for a column that == valueAsCol
        function idx = findCol(obj, valueAsCol)
            
            idx = find(obj.values(1,:) == valueAsCol(1));
            
            if numel(idx) > 1
                dims = numel(valueAsCol);
                for dd = 1:dims
                    idx = idx(obj.values(dd,idx) == valueAsCol(dd));
                    
                    if numel(idx) == 1
                        return;
                    end
                end
            end
            
            if numel(idx) > 1
                warning('Found multiple.');
                idx = idx(1);
            end
            
        end

        function idx = locate(obj, value)
            if obj.isRow
                assert(iscolumn(value));
                idx = obj.findCol(value);
            else
                assert(isrow(value));
                idx = obj.findRow(value);
            end
        end
        
        
        % For the heap element with the given value,
        % replace its key and re-heapify.
        function reKey(obj, value, key)
            
            idx = obj.locate(value);
            obj.reKeyIndex(idx, key);
            
        end
        
        % Given the particular index, re-key it.
        function reKeyIndex(obj, idx, key)
            obj.clearSwapBuffer();
            
            if idx ~= 0
                obj.keys(idx) = key;
                
                parent = @(n) floor(n/2);
                
                if idx > 1 && key > obj.keys(parent(idx))
                    obj.bubbleUp(idx);
                else
                    obj.bubbleDown(idx);
                end
            end
        end
            
        
    end % methods
    
    
end