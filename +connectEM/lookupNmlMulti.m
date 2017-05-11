function [segIds, neighbours, filenames, nodes, startNode, comments] = lookupNmlMulti(p, folders, removeNmlWithMoreThanOneTree)

    if nargin < 3
        removeNmlWithMoreThanOneTree = true;
    end


    tic;
    for f=1:length(folders)
        files = dir([folders{f} '*.nml']);
        filenames{f} = cellfun(@(x)[folders{f} x], {files(:).name}, 'uni', 0);
    end
    toc;
    tic;
    filenames = cat(2, filenames{:});
    idxTooManyTrees = false(length(filenames),1);
    comments = cell(length(filenames),1);
    for i=1:length(filenames)
        try
            [~,skel] = evalc('parseNml_webKnossos(filenames{i})');
        catch ME
            display(filenames{i});
            warning('parseNml generated an error');
            display(ME.message);
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
            nodes{i} = skel{1}.nodesNumDataAll(:,3:5);
            time = skel{1}.nodesNumDataAll(:,8);
            [~, idxT] = sort(time);
            startNode{i} = nodes{i}(idxT(1),:);
            % Fix nodes (double nodes at same position)
            nodes{i} = unique(nodes{i}, 'rows');
        end
    end
    filenames(idxTooManyTrees) = [];
    nodes(idxTooManyTrees) = [];
    startNode(idxTooManyTrees) = [];
    toc;
    temp.nodes = nodes;
    [segIds, neighbours] = connectEM.lookupSkelGT(p, temp);

end

