function structToHdf5(file, group, struct, forceArray)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    assert(isstruct(struct));
    
    if ~exist('forceArray', 'var') || isempty(forceArray)
        forceArray = false;
    end
    
    if ~isscalar(struct) || forceArray
        for curIdx = 1:numel(struct)
            curGroup = fullfile(group, num2str(curIdx));
            structToHdf5(file, curGroup, struct(curIdx));
        end
    else
        names = fieldnames(struct);
        values = struct2cell(struct);
        assert(isequal(numel(names), numel(values)));

        for curIdx = 1:numel(names)
            curDset = fullfile(group, names{curIdx});
            arrayToHdf5(file, curDset, values{curIdx});
        end
    end
end
