function jobs = globalization(p)

    display('Constructing array for segmentation ID offset in cubes');
    tic;
    % collect number of segments per cube
    cumNum = uint32(0);
    numElTotal = zeros(size(p.local), 'uint32'); 
    numElTotalUpper = zeros(size(p.local), 'uint32');
    for i=1:size(p.local,1)
        for j=1:size(p.local,2)
            for k=1:size(p.local,3)
                load(p.local(i,j,k).segFile, 'numEl');
                numElTotal(i,j,k) = cumNum;
                cumNum = cumNum + uint32(numEl);
                numElTotalUpper(i,j,k) = numEl+1; 
            end
        end
    end
    % Save numElTotal so that it only has to be added to localID of respective cube to get global one
    Util.save([p.saveFolder 'numEl.mat'], numElTotal, numElTotalUpper); 
    toc;

    % Globalize segmentation and save as KNOSSOS hierachy for Oxalis
    for i=1:size(p.local,1)
        for j=1:size(p.local,2)
            for k=1:size(p.local,3)
                idx = sub2ind(size(p.local), i, j, k);
                inputCell{idx} = {p, i, j, k};
            end
        end
    end
    functionH = @globalSegId;
    jobs(1) = startCPU(functionH, inputCell, 'globalSegmentID', 12, 10);

    % Globalize correspondences
    clear functionH inputCell;
    files = dir([p.correspondence.saveFolder, '*.mat']);
    % Bad but internal: Sort out global correspondences in case pipeline has already run
    idx = cellfun(@length, {files(:).name}) == 16;
    files = files(idx);
    for i=1:length(files)
        inputCell{i} = {p, files(i).name};
    end
    functionH = @globalCorrSeg;
    jobs(2) = startCPU(functionH, inputCell, 'globalCorrespondences', 12, 100);

end
