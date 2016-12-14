function job = globalizeSegmentation(p)
    display('Collecting segment ID offsets...');

    % collect number of segments per cube
    cumNum = uint32(0);
    numElTotal = zeros(size(p.local), 'uint32'); 
    numElTotalUpper = zeros(size(p.local), 'uint32');

    tic();
    for curIdx = 1:numel(p.local)
        load(p.local(curIdx).segFile, 'numEl');

        numElTotal(curIdx) = cumNum;
        cumNum = cumNum + uint32(numEl);
        numElTotalUpper(curIdx) = numEl + 1;
    end
    toc();

    % Save numElTotal so that it only has to be added to localID of respective cube to get global one
    Util.save([p.saveFolder 'numEl.mat'], numElTotal, numElTotalUpper); 

    % collect parameters
    inputCell = cell(numel(p.local), 1);

    for curIdx = 1:numel(p.local)
        [curI, curJ, curK] = ind2sub( ...
            size(p.local), curIdx);
        inputCell{curIdx} = {p, curI, curJ, curK};
    end

    % Globalize segmentation and save as KNOSSOS hierachy for Oxalis
    job = startCPU(@globalSegId, inputCell, 'globalSegmentID', 12, 10);
end

