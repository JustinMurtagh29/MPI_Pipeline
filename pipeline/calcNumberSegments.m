function calcNumberSegments(p)
    % Calculate offset local to global segmentation for each cube

    coords = localToGlobalBboxModification(p);

    numEl = uint32(0);
    numElTotal = zeros(size(p.local), 'uint32');
    for i=1:size(p.local,1)
        for j=1:size(p.local,2)
            for k=1:size(p.local,3)
                % Load segmentation and extract relevant section
                load(p.local(i,j,k).segFile)
                seg = seg(coords{i,j,k}(1,1):coords{i,j,k}(1,2), coords{i,j,k}(2,1):coords{i,j,k}(2,2), coords{i,j,k}(3,1):coords{i,j,k}(3,2));
                % How many unique global IDs are needed for this cube
                globalIds = unique(seg(:));
                nrGlobalIdsNeededForThisCube = length(globalIds);
                % If there is a 0, do not count this one 
                if any(globalIds == 0)
                    nrGlobalIdsNeededForThisCube = nrGlobalIdsNeededForThisCube - 1;
                end
                % Update sum of all former nrGlobalIDs & collect for all cubes for use with correspondences
                numElTotal(i,j,k) = numEl;
                numEl = numEl + uint32(nrGlobalIdsNeededForThisCube);
            end
        end
    end
    % Save numElTotal so that it only has to be added to localID of respective cube to get global one
    save([p.seg.root 'numEl.mat'], 'numElTotal'); 

end

