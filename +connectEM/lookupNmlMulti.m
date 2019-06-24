function [segIds, neighbours, filenames, nodes, startNode, comments, errors] = lookupNmlMulti(p, folders, removeNmlWithMoreThanOneTree)
% loads all nmls within the locations descriped in "folders"
    if nargin < 3
        removeNmlWithMoreThanOneTree = true;
    end

    tic;
    filenames = cell(1, numel(folders));
    for f=1:length(folders)
        files = dir(fullfile(folders{f}, '*.nml'));
        filenames{f} = cellfun(@(x) fullfile(folders{f}, x), {files(:).name}, 'uni', 0);
    end
    toc;
    
    tic;
    filenames = cat(2, filenames{:});
    nodes = cell(1, numel(filenames));
    startNode = cell(1, numel(filenames));
    idxTooManyTrees = false(length(filenames),1);
    comments = cell(length(filenames),1);
    errors = cell(length(filenames),1);
    
    for i=1:length(filenames)
        try
            [~,skel] = evalc('parseNml_webKnossos(filenames{i})');
        catch ME
            display(filenames{i});
            warning('parseNml generated an error');
            display(ME.message);
            errors{i,1}{1} = ME.message;
            errors{i,1}{2} = filenames{i};
            continue;
        end
        
        if length(skel) >= 2
            if removeNmlWithMoreThanOneTree
                display(filenames{i});
                warning('Skeleton has too many trees');
                idxTooManyTrees(i) = true;
            else
                nodes{i} = [];
                for j=1:length(skel)
                    nodes{i}(end+1:end+size(skel{j}.nodesNumDataAll,1),1:3) = skel{j}.nodesNumDataAll(:,3:5);
                end
                nodes{i} = unique(nodes{i}, 'rows');
                % Start node not yet implemented if multiple trees in nml (needed for queries, but not GT)
                startNode{i} = [];
            end
        else
            if ~isempty(skel{1}.commentsString)
                comments{i} = skel{1}.commentsString;
            end
            % Extract needed information from skeletons
			disp (['processing file: ' num2str(filenames{i})])
            nodes{i} = skel{1}.nodesNumDataAll(:,3:5);
            time = skel{1}.nodesNumDataAll(:,8);
            [~, idxT] = sort(time);
            startNode{i} = nodes{i}(idxT(1),:);
            nodes{i} = nodes{i}(idxT,:);
            % Fix nodes (double nodes at same position)
            nodes{i} = unique(nodes{i}, 'rows','stable'); % preserve the order of time
        end
    end
    
    filenames(idxTooManyTrees) = [];
    nodes(idxTooManyTrees) = [];
    startNode(idxTooManyTrees) = [];
    comments(idxTooManyTrees) = [];
    errors(idxTooManyTrees) = [];
    toc;
    
    if ~isempty(p)
        temp.nodes = nodes;
        [segIds, neighbours] = connectEM.lookupSkelGT(p, temp);
    else
        segIds = [];
        neighbours = [];
    end

end
