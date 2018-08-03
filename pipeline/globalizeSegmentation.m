function job = globalizeSegmentation(p)
    Util.log('Collecting segment ID offsets...');

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

    % Save numElTotal so that it only has to be added to localID of
    % respective cube to get global one
    Util.save([p.saveFolder 'numEl.mat'], numElTotal, numElTotalUpper); 

    % collect parameters
    inputCell = cell(numel(p.local), 1);
    for curIdx = 1:numel(p.local)
       [curI, curJ, curK] = ind2sub(size(p.local), curIdx);
        inputCell{curIdx} = {curI, curJ, curK};
    end
    
    % init wkw dataset, if needed
    if isfield(p.seg, 'backend') && strcmp(p.seg.backend,'wkwrap')
        wkwInit('new', p.seg.root, 32, 32, 'uint32', 1);
    end
    
    job = Cluster.startJob( ...
        @globalSegId, inputCell, ...
        'sharedInputs', {p}, ...
        'name', mfilename(), ...
        'taskGroupSize', 10, ...
        'cluster', {'memory', 12, 'time','04:00:00'});
end
