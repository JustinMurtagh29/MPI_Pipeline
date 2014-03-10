function skel = removeEmptySkeletons( skel )

toDel = zeros(size(skel,2),1);
for i=1:length(skel)
    if isempty(skel{i}.nodes) || size(skel{i}.nodes,1) == 0
        toDel(i) = 1;
    end
end
skel = skel(~toDel);
display(['Removed ' num2str(sum(toDel)) ' empty skeletons.']);

end

