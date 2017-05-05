function [er, cm] = getERcomponents()
% Eine Perle der Codegeschichte

    % first round of ER annotation
    temp = load('+connectEM/evaluationData/gliaER/00022.mat', 'er', 'axonsFinal');
    gliaERannotations = skeleton('+connectEM/evaluationData/gliaER/2012-09-28_ex145_07x2_ROI2017__explorational__mberning__a361bb_gliaER.nml');
    treeNameInfo = cellfun(@(x)regexp(x, 'axonsLargest(\d\d\d\d)_.*_15030572( - gliaER)', 'tokens'), gliaERannotations.names, 'uni', 0);
    gliaERidx = cell2mat(cellfun(@(x)str2num(x{1}{1}), treeNameInfo(~cellfun(@isempty, treeNameInfo)), 'uni', 0));
    % second round of ER (and CM) annotation
    temp2 = load('+connectEM/evaluationData/gliaER/cmAndErAnnotations.mat');
    gliaERannotations2 = skeleton('+connectEM/evaluationData/gliaER/2012-09-28_ex145_07x2_ROI2017__explorational__mberning__d85b75.nml');
    treeNameInfo2 = cellfun(@(x)regexp(x, 'axonsLargest93_(\d\d\d\d)_.*_15030572( - remainingER)', 'tokens'), gliaERannotations2.names, 'uni', 0);
    gliaERidx2 = cell2mat(cellfun(@(x)str2num(x{1}{1}), treeNameInfo2(~cellfun(@isempty, treeNameInfo2)), 'uni', 0));
    er = cat(1, temp.er, temp.axonsFinal(gliaERidx), temp2.axonsFinal(gliaERidx2));
    % special section for catastrpohic merger extraction
    treeNameInfo3 = cellfun(@(x)regexp(x, 'axonsLargest93_(\d\d\d\d)_.*_15030572( - cm)', 'tokens'), gliaERannotations2.names, 'uni', 0);
    gliaERidx3 = cell2mat(cellfun(@(x)str2num(x{1}{1}), treeNameInfo3(~cellfun(@isempty, treeNameInfo3)), 'uni', 0));
    cm = temp2.axonsFinal(gliaERidx3);

end

