function [segIds, neighbours, filenames, nodes] = lookupNml(p, folder)

    files = dir([folder '*.nml']);
    filenames = {files(:).name};
    for i=1:length(files)
        [~,skel{i}] = evalc('parseNml_webKnossos([folder files(i).name])');
        nodes{i} = skel{i}{1}.nodes(:,1:3);
        [~, idx] = unique(nodes{i}, 'rows');
        nodes{i} = nodes{i}(idx,:);
    end
    temp.nodes = nodes;
    [segIds, neighbours] = lookupSkelGT(p, temp);

end

