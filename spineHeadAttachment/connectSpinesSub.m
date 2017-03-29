function good_list = connectSpinesSub(starting, stopIDs, ...
    connM, maxSteps, exclusions)
for idx = 1 : length(starting)
    idx
    good_list.comment{idx} = {makeComment('start', exclusions, starting(idx))};
    good_list.segIds{idx} = starting(idx);
    good_list.connectTo{idx} = -1;
    while isempty(intersect(stopIDs, good_list.segIds{idx})) ...
            && length(good_list.segIds{idx}) < maxSteps
        neighbours = connM(:, good_list.segIds{idx});
        neighbours(good_list.segIds{idx}, :) = 0;
        [value, idx2] = max(neighbours(:));
        if ~value
            break;
        end
        [newid, oldid] = ind2sub(size(neighbours),idx2);
        
        good_list.comment{idx}{end + 1} = ...
            makeComment(num2str(connM(good_list.segIds{idx}(oldid), newid)), exclusions, newid);
        good_list.segIds{idx} = [good_list.segIds{idx} newid];
        good_list.connectTo{idx}(end + 1) = find(good_list.segIds{idx} == good_list.segIds{idx}(oldid));
        fprintf('step %u\n', newid);
    end
end
