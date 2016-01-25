function skelGT = skeletons07x2()
    skelFolder = '/gaba/u/mberning/2012-09-28_ex145_07x2_skeletons/';

    % Whole Cell Tracings
    wholeCellSubfolder = [skelFolder 'somaSeeds/'];
    files = dir([wholeCellSubfolder 'MINK*']);
    for i=1:length(files)
        files2 = dir([wholeCellSubfolder files(i).name '/all/' files(i).name '/*.nml']);
        for j=1:length(files2)
            fieldname = ['id' num2str(i, '%.4i')]; 
            fieldname2 = ['skel' num2str(j, '%.4i')];
            skelGT.wholeCell.(fieldname).(fieldname2).file = [wholeCellSubfolder files(i).name '/all/' files(i).name '/' files2(j).name];
            skelGT.wholeCell.(fieldname).(fieldname2).skel = skeleton(skelGT.wholeCell.(fieldname).(fieldname2).file, false, 1);
        end
    end

    % Apicals
    skelGT.apicals.file = [skelFolder 'apicals/' 'apikaldendritenfinal25.012_resaved_firstcleanup_postmatlab.nml'];
    skelGT.apicals.skel = skeleton(skelGT.apicals.file, false, 1);

    % Dense Dendrite Tracings
    denseDendriteSubfolder = [skelFolder 'dendrites/'];
    files = dir([denseDendriteSubfolder '*.nml']);
    for i=1:length(files)
        fieldname = ['id' num2str(i, '%.4i')];
        skelGT.denseDendrite.(fieldname).file = [denseDendriteSubfolder files(i).name];
        skelGT.denseDendrite.(fieldname).skel = skeleton(skelGT.denseDendrite.(fieldname).file, false, 1);
    end

    % Inital Segment Innervation Tracings
    skelGT.axonsIS.file = [skelFolder 'axonsFromIS/' 'output_final_all_joint_oxalisoutput_iwanttodiffthis.nml'];
    skelGT.axonsIS.skel = skeleton(skelGT.axonsIS.file, false, 1);

    % Dendrite Innervation Tracings
    dendriteInnervationSubfolder = [skelFolder 'post_iris/'];
    files = dir([dendriteInnervationSubfolder '*.nml']);
    for i=1:length(files)
        fieldname = ['axonsFromDI' num2str(i, '%.4i')];
        skelGT.(fieldname).file = [dendriteInnervationSubfolder files(i).name];
        skelGT.(fieldname).skel = skeleton(skelGT.(fieldname).file, false, 1);
        % Remove 3-Point annotation from SynTypComm file
        if ~isempty(strfind(skelGT.(fieldname).file, 'SynTypComm'));
            for i=1:length(skelGT.(fieldname).skel)
                skelGT.(fieldname).skel = removeTPAsSimple(skelGT.(fieldname).skel,i);
            end
        end
    end

    % Oblique Spine Innervation Tracings
    % Dense Dendrite Tracings
    osiSubfolder = [skelFolder 'obliqueSpineInnervation/'];
    files = dir([osiSubfolder '*.nml']);
    for i=1:length(files)
        fieldname = ['id' num2str(i, '%.4i')];
        skelGT.axonsFromOS.(fieldname).file = [osiSubfolder files(i).name];
        skelGT.axonsFromOS.(fieldname).skel = skeleton(skelGT.axonsFromOS.(fieldname).file, false, 1);
    end

end

