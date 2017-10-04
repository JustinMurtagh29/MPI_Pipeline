function [cc, centerOfCC] = ...
        detectChiasmataNodesCluster(nodes, isIntersection)
    if sum(isIntersection) > 1
        N = sum(isIntersection);
        if N > 10000
            distances = pdist(nodes(isIntersection,:));
            connections = [];
            for idx = 1 : 100
                connections = [connections, find(distances(idx:100:end)<2000)];
            end
            clear distances
            connections = sort(connections);
            % pdist lists the distances in a special way (see documentation)
            % the next two lines define the borders of the right column of those replies
            idxs_end = cumsum(N-1:-1:1);
            idxs_start = [1, cumsum(N-1:-1:2)+1];
            % sort the connections within these boundaries
            [B, I] = sort([connections,(idxs_end + 0.5)]); % B = A(I);
            I_reverse(I) = 1:length(I);
            % find the separators again in the sorted list
            separators = find(B ~= round(B));
            % make a lookup to find for each entry in which box it is
            bins = repelem(1:(N-1),diff([1,separators]));
            %with that we get the right side for each connection
            rightside = bins(I_reverse(1 : end-1));
            rightside = rightside(1:length(connections));
            %the left side depends where ipwdt is in its box
            leftside_pre = idxs_start(rightside);
            leftside = connections-leftside_pre+rightside+1;
            [cccount,ccpre] =graphconncomp(sparse(leftside, rightside,ones(1,length(connections)),N,N),'Directed',false);
        else
            % NOTE(amotta): Group chiasmatic nodes with less than 2 µm distance
            [cccount,ccpre] =graphconncomp(sparse(squareform(pdist(nodes(isIntersection,:)))<2000));
        end
        lookup = find(isIntersection);
        cc = arrayfun(@(x){lookup(ccpre==x)},1:cccount);

        % NOTE(amotta): Pick node closest to c.o.m. as center
        [~, centerOfCC] = cellfun(@(x)min(pdist2(bsxfun(@minus, nodes(x,:), mean(nodes(x,:),1)), [0 0 0])), cc);
    else 
        if sum(isIntersection) == 1
            cc = {find(isIntersection)};
            centerOfCC = 1;
        else
            cc = {};
            centerOfCC = [];
        end
    end
end