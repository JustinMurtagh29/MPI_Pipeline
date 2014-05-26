function leaves = findLeaves(edges, segIDsToBorder)
%In leaves{k,1} steht die ObjectId der "Muttersupervoxel", in leaves{k,2} die
%umschlossenen Supervoxel

objectIds = unique(edges);
k = 1;
leaves = {};

for i=1:size(objectIds,1)
    
    %Herausloeschen des Objekts i und Graph auf Zusammenhangskomponenten
    %testen
    newedgelist = edges;
    index = any(newedgelist == objectIds(i), 2);
    newedgelist (index,:)= [];
    newobjectIds = objectIds;
    newobjectIds(i) = [];
    components = findConnectedComponents(newedgelist,newobjectIds);
    
    %Wenn mehr als 1 Zusammenhangskomponente: Verdacht auf leaves, objekte
    %sind aber nur mit Sicherheit leaves, wenn alle Objekte am border
    %liegen
    if size(components,1)>1
        for j= 1:size(components,1)
            if isempty(intersect(components{j},segIDsToBorder))
                leaves{k,1} = objectIds(i);
                leaves{k,2} = components{j};
                k = k+1;
            end
        end
    end
end

end
