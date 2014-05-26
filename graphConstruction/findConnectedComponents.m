function components = findConnectedComponents(newedgelist,newobjectIds)

components ={};
k = 1;
%Liste der noch nicht besuchten Knoten
nonvisited = 1:size(newobjectIds,1);

while ~isempty(nonvisited)
    components{k,1} = newobjectIds(nonvisited(1));
    l1 = size(components{k,1},1);
    %Hinzufuegen aller Nachbarn von nonvisited(1)
    index = any(newedgelist == newobjectIds(nonvisited(1)), 2);
    v = newedgelist(index,:);
    %Loeschen des besuchten Wertes
    nonvisited(1)=[];
    v = intersect(v(:),newobjectIds(nonvisited));
    components{k,1}=vertcat(components{k,1},v);
    l2 = size(components{k,1},1);
    %Hinzufuegen aller Nachbarn der Nachbarn etc
    while ~(l1==l2)
        v=[];
        for j=l1+1:l2
            index = find(newobjectIds==(components{k,1}(j)));
            index1 = any(newedgelist == newobjectIds(index), 2);
            v=[v; newedgelist(index1,:)];
            nonvisited(nonvisited==index)=[];
        end
        l1 = l2;
        v = intersect(v(:),newobjectIds(nonvisited));
        components{k,1}=vertcat(components{k,1},v);
        l2 = size(components{k,1},1);
    end
    k = k+1;
end

end