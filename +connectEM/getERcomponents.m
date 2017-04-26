function er = getERcomponents()
% Eine Perle der Codegeschichte

    temp = load('+connectEM/evaluationData/gliaER/00022.mat', 'er', 'axonsFinal');
    gliaERannotations = skeleton('+connectEM/evaluationData/gliaER/2012-09-28_ex145_07x2_ROI2017__explorational__mberning__a361bb_gliaER.nml');
    treeNameInfo = cellfun(@(x)regexp(x, 'axonsLargest(\d\d\d\d)_.*_15030572( - gliaER)', 'tokens'), gliaERannotations.names, 'uni', 0);
    gliaERidx = cell2mat(cellfun(@(x)str2num(x{1}{1}), treeNameInfo(~cellfun(@isempty, treeNameInfo)), 'uni', 0));
    er = cat(1, temp.er, temp.axonsFinal(gliaERidx));

end

