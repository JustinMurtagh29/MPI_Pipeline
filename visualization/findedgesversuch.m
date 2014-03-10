function [edges, borders] = findedgesversuch(seg, conn)

maxId = max(seg(:))+1;
segpad = padarray(seg,[1 1 1],maxId);
[M,N,P] = size(segpad);

if conn == 6
    vec = [-1 1 -M M -M*N +M*N];
elseif conn == 26
    vec = [(-M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1]) [-M-1 -M -M+1 -1 1 M-1 M M+1] (M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1])];
end

ind = find(segpad==0);
nInd = bsxfun(@plus, ind, vec);
nSegId = reshape(segpad(nInd(:)), size(nInd));
nSegCell = mat2cell(nSegId, ones(size(nSegId,1),1), conn);
fun1A = @(x)sort(x);
fun1B = @(x)x([true diff(x) > 0]);
list1 = cellfun(@(x)nonzeros(fun1B(fun1A(x))),nSegCell,'UniformOutput', false);
fun2 = @(x) combnk(x,2);
pairs = cellfun(fun2,list1,'UniformOutput', false);
pairs = cell2mat(pairs);
edges = unique(pairs, 'rows');
edges(edges(:,2) == maxId,:) = [];

end
