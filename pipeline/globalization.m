function globalization(parameter)
    % calculates number of segments per cube

    coords = localToGlobalBboxModification(parameter);
   
    % Takes approx 3 hours on 07x2, paralellize? 
    display('Calculating ID offset local to global bounding box');
    tic;
    calcNumberSegments(parameter, coords);
    toc;

    % Globalize segmentation and save as KNOSSOS hierachy for Oxalis
    for i=1:size(parameter.local,1)
        for j=1:size(parameter.local,2)
            for k=1:size(parameter.local,3)
                idx = sub2ind(size(parameter.local), i, j, k);
                functionH{idx} = @globalSegId;
                inputCell{idx} = {parameter, coords, i, j, k};
            end
        end
    end
    startCPU(functionH, inputCell, 'globalize segmentation IDs');

    % Globalize correspondences
    clear functionH inputCell;
    files = dir([parameter.correspondence.saveFolder, '*.mat']);
    for i=1:length(files)
        functionH{i} = @globalCorrSeg;
        inputCell{i} = {parameter, files(i).name};
    end
    startCPU(functionH, inputCell, 'globalize correspondences');

end
