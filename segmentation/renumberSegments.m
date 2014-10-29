function stackBig = renumberSegments( stackBigWithoutRenumbering, stackArray, threshold )
            
planeOld = cell(numel(stackArray),6);
planeNew = cell(numel(stackArray),6);
for x=1:size(stackArray,1)
    for y=1:size(stackArray,2)
        for z=1:size(stackArray,3)
            lowerBound = ([x y z] - 1)*100+1;
            upperBound = [x y z]*100;
            if lowerBound(1) > 1
                planeOld{sub2ind(size(stackArray),x,y,z),1} = squeeze(stackBigWithoutRenumbering(lowerBound(1)-1,lowerBound(2):upperBound(2),lowerBound(3):upperBound(3)));
                planeNew{sub2ind(size(stackArray),x,y,z),1} = squeeze(stackArray{x,y,z}(1,:,:));  
            end
            if lowerBound(2) > 1
                planeOld{sub2ind(size(stackArray),x,y,z),2} = squeeze(stackBigWithoutRenumbering(lowerBound(1):upperBound(1),lowerBound(2)-1,lowerBound(3):upperBound(3)));
                planeNew{sub2ind(size(stackArray),x,y,z),2} = squeeze(stackArray{x,y,z}(:,1,:));
            end
            if lowerBound(3) > 1
                planeOld{sub2ind(size(stackArray),x,y,z),3} = squeeze(stackBigWithoutRenumbering(lowerBound(1):upperBound(1),lowerBound(2):upperBound(2),lowerBound(3)-1));
                planeNew{sub2ind(size(stackArray),x,y,z),3} = squeeze(stackArray{x,y,z}(:,:,1));
            end
            if upperBound(1) < size(stackBigWithoutRenumbering,1)
                planeOld{sub2ind(size(stackArray),x,y,z),4} = squeeze(stackBigWithoutRenumbering(upperBound(1)+1,lowerBound(2):upperBound(2),lowerBound(3):upperBound(3)));
                planeNew{sub2ind(size(stackArray),x,y,z),4} = squeeze(stackArray{x,y,z}(end,:,:));
            end
            if upperBound(2) < size(stackBigWithoutRenumbering,2)
                planeOld{sub2ind(size(stackArray),x,y,z),5} = squeeze(stackBigWithoutRenumbering(lowerBound(1):upperBound(1),upperBound(2)+1,lowerBound(3):upperBound(3)));
                planeNew{sub2ind(size(stackArray),x,y,z),5} = squeeze(stackArray{x,y,z}(:,end,:));
            end
            if upperBound(3) < size(stackBigWithoutRenumbering,3)
                planeOld{sub2ind(size(stackArray),x,y,z),6} = squeeze(stackBigWithoutRenumbering(lowerBound(1):upperBound(1),lowerBound(2):upperBound(2),upperBound(3)+1));
                planeNew{sub2ind(size(stackArray),x,y,z),6} = squeeze(stackArray{x,y,z}(:,:,end));
            end
        end
    end
end
planeNew(cellfun(@isempty,planeNew)) = [];
planeOld(cellfun(@isempty,planeOld)) = [];

uniqueIDsOld = unique(stackBigWithoutRenumbering);
uniqueIDsOld(uniqueIDsOld == 0) = [];
joinIDs = zeros(length(uniqueIDsOld));
for plane=1:length(planeOld)
    for i=1:length(uniqueIDsOld)
        for j=1:i
            overlap = (planeNew{plane} == uniqueIDsOld(i)) & (planeOld{plane} == uniqueIDsOld(j));
            if sum(overlap(:)) > threshold
                joinIDs(i,j) = 1;
            end
        end
    end
end

stackBig = stackBigWithoutRenumbering;
for i=1:length(uniqueIDsOld)
    queue = i;
    newLabels = [];
    while ~isempty(queue)
        newLabels = [newLabels queue(1)];
        queue = unique([queue find(joinIDs(queue(1),:)) find(joinIDs(:, queue(1))')]);
        for j=1:length(newLabels)
            queue(queue == newLabels(j)) = [];
        end
    end
    for j=1:length(newLabels)
        stackBig(stackBig == uniqueIDsOld(newLabels(j))) = uniqueIDsOld(i);
    end
end

end
