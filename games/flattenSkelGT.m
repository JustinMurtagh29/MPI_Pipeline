function skelGTflat = flattenSkelGT(skelGT)

    % Flatten structure into cell arrays
    fNames = fieldnames(skelGT);
    skelGTflat.filename = cell(0);
    skelGTflat.treename = cell(0);
    skelGTflat.nodes = cell(0);
    skelGTflat.edges = cell(0);
    idx = 1;  
    for i=1:length(fNames)
        % Flatten structure, different kind of input formats, here files in deep folder hierachy
        if strcmp(fNames{i}, 'wholeCell')
            fNames2 = fieldnames(skelGT.(fNames{i}));
            for j=1:length(fNames2)
                fNames3 = fieldnames(skelGT.(fNames{i}).(fNames2{j}));
                for k=1:length(fNames3)
                    for l=1:length(skelGT.(fNames{i}).(fNames2{j}).(fNames3{k}).skel)
                        skelGTflat.filename{idx} = skelGT.(fNames{i}).(fNames2{j}).(fNames3{k}).file;
                        skelGTflat.treename{idx} = skelGT.(fNames{i}).(fNames2{j}).(fNames3{k}).skel.names{l};
                        skelGTflat.nodes{idx} = skelGT.(fNames{i}).(fNames2{j}).(fNames3{k}).skel.nodes{l}(:,1:3);
                        skelGTflat.edges{idx} = skelGT.(fNames{i}).(fNames2{j}).(fNames3{k}).skel.edges{l};
                        idx = idx + 1;
                    end
                end
            end
        % Here all nml in one folder
        elseif strcmp(fNames{i}, 'denseDendrite') || strcmp(fNames{i}, 'axonsFromOS');
            fNames2 = fieldnames(skelGT.(fNames{i}));
            for j=1:length(fNames2)
                for k=1:length(skelGT.(fNames{i}).(fNames2{j}).skel)
                    skelGTflat.filename{idx} = skelGT.(fNames{i}).(fNames2{j}).file;
                    skelGTflat.treename{idx} = skelGT.(fNames{i}).(fNames2{j}).skel.names{k};
                    skelGTflat.nodes{idx} = skelGT.(fNames{i}).(fNames2{j}).skel.nodes{k}(:,1:3);
                    skelGTflat.edges{idx} = skelGT.(fNames{i}).(fNames2{j}).skel.edges{k};
                    idx = idx + 1;
                end
            end
        % Here for different trees in one file
        else
            for j=1:length(skelGT.(fNames{i}).skel.names)
                skelGTflat.filename{idx} = skelGT.(fNames{i}).file;
                skelGTflat.treename{idx} = skelGT.(fNames{i}).skel.names{j};
                skelGTflat.nodes{idx} = skelGT.(fNames{i}).skel.nodes{j}(:,1:3);
                skelGTflat.edges{idx} = skelGT.(fNames{i}).skel.edges{j};
                idx = idx + 1;
            end
        end
    end

end

