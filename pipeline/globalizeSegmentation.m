
function job = globalizeSegmentation(p)
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
    save([p.saveFolder 'numEl.mat'], 'numElTotal', 'numElTotalUpper'); 
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
    job = startCPU(functionH, inputCell, 'globalSegmentID', 12, 10);

end

