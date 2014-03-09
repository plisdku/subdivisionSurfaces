% Test the heap!

import VVMesh.*

%% Build a heap from scratch.

heap = Heap;

heap.push(1);
heap.push(3);
heap.push(2);
heap.push(5);
heap.push(4);
heap.push(6);
assert(heap.heapSize == 6);

popVals = sort(heap.keys, 'descend'); % Expected order to pop values in

heapSize = numel(heap);
for ii = 1:numel(popVals)
    assert(heap.pop() == popVals(ii));
end

assert(heap.heapSize == 0);

fprintf('PASSED pop test with %i pops!\n', numel(popVals));

%% Build a heap from an array.

myArray = rand(10,1);

heap = Heap(myArray);
assert(heap.heapSize == numel(myArray));

sorted = sort(myArray, 'descend');

for ii = 1:numel(sorted)
    assert(heap.pop() == sorted(ii));
end
assert(heap.heapSize == 0);

fprintf('PASSED constructor test!\n');

%% Data test with a string

heap = Heap;

myString = 'commando';
myKeys = numel(myString):-1:1;
p = randperm(numel(myString));

for ii = 1:numel(myString)
    heap.push(myKeys(p(ii)), myString(p(ii)));
end

fprintf('Scrambled "%s" to "%s" and inserted in heap.\n',...
    myString, myString(p));

reStr = '';

for ii = 1:numel(myString)
    [~, reStr(ii)] = heap.pop();
end

assert(strcmp(myString, reStr));

fprintf('PASSED string test!\n');

%% Data test with arrays

strings = ['is'; 'it'; 'my'; 'ka'];
scores = [1; 5; 3; 7];

heap = Heap(scores, strings);

while heap.heapSize > 0
    [~, str] = heap.pop();
    fprintf('%s\n', str);
end