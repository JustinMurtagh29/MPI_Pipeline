function cc = findCC(edges)

% To ignore mutliple borders between objects in construction of connectivity matrix
edges = unique(sort(edges,2), 'rows');
% to get a symetric sparse matrix
edges = [edges; edges(:,[2 1])];

% Using Dual Mendelsohn decomposition for finding connected components 
n = max(edges(:));
A = sparse([edges(:,1); (1:n)'],[edges(:,2); (1:n)'],1);
[p,q,r,s] = dmperm(A);
cc = cell(length(r)-1,1);
for i=1:length(r)-1
	if (r(i+1)-r(i)) > 1
		cc{i} = p(r(i):r(i+1)-1);
	end
end
cc(cellfun(@isempty,cc)) = [];

end

