function label = findLabels( edges, trueEdges )

label = zeros(size(edges,1),1);
for i=1:size(edges,1)
    if any(sum(repmat(edges(i,:),[size(trueEdges,1) 1]) == trueEdges,2) == 2)
        label(i) = 1;
    else
	label(i) = -1;
    end
end

end

