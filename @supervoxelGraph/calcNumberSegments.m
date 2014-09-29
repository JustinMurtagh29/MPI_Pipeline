function calcNumberSegments(p)

% construct coord array for bBox border area for each cube in local field of p
coords{1} = [257,768;257,768;129,384];
coords = repmat(coords,[size(p.local,1),size(p.local,2),size(p.local,3)]);
% Usually just use inner bounding box, at outer limits of ROI borders needed to always have enough segmentation context
coords(1,:,:) = cellfun(@(x)([x(1,1)+p.tileBorder(1,1) x(1,2); x(2,1) x(2,2); x(3,1) x(3,2)]),coords(1,:,:), 'UniformOutput', false);
coords(end,:,:) = cellfun(@(x)([x(1,1) x(1,2)+p.tileBorder(1,2); x(2,1) x(2,2); x(3,1) x(3,2)]),coords(end,:,:), 'UniformOutput', false);
coords(:,1,:) = cellfun(@(x)([x(1,1) x(1,2); x(2,1)+p.tileBorder(2,1) x(2,2); x(3,1) x(3,2)]),coords(:,1,:), 'UniformOutput', false);
coords(:,end,:) = cellfun(@(x)([x(1,1) x(1,2); x(2,1) x(2,2)+p.tileBorder(2,2); x(3,1) x(3,2)]),coords(:,end,:), 'UniformOutput', false);
coords(:,:,1) = cellfun(@(x)([x(1,1) x(1,2); x(2,1) x(2,2); x(3,1)+p.tileBorder(3,1) x(3,2)]),coords(:,:,1), 'UniformOutput', false);
coords(:,:,end) = cellfun(@(x)([x(1,1) x(1,2); x(2,1) x(2,2); x(3,1) x(3,2)+p.tileBorder(3,2)]),coords(:,:,end), 'UniformOutput', false);

% globalize segmentaiton Ids
numEl = uint32(0);
numElTotal = zeros(size(p.local), 'uint32');
numElTotalAll = zeros(size(p.local), 'uint32');
for i=1:size(p.local,1)
    for j=1:size(p.local,2)
        for k=1:size(p.local,3)
            % Load segmentation and extract relevant section
            load(p.local(i,j,k).segFile)
            seg = seg(coords{i,j,k}(1,1):coords{i,j,k}(1,2), coords{i,j,k}(2,1):coords{i,j,k}(2,2), coords{i,j,k}(3,1):coords{i,j,k}(3,2));
            % Remember what is 0 and how many global IDs are needed for this cube
            isZero = seg == 0;
            nrGlobalIDs = length(unique(seg(:)))-1;
            % Update sum of all former nrGlobalIDs & collect for all cubes for use with correspondences
            numElTotal(i,j,k) = numEl;
            numEl = numEl + uint32(nrGlobalIDs);
            numElTotalAll(i,j,k) = numEl;
        end
    end
end
% Save numElTotal so that it only has to be added to localID of repective cube to get global one
save([p.seg.root 'numEl.mat'], 'numElTotal') 

end

