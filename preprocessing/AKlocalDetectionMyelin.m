function AKlocalDetectionMyelin(in, out, bboxIn, bboxOut)

    % Load raw data
    raw = readKnossosRoi(in.root, in.prefix, bboxIn);
    % Load blood vessel segmentation
    
    vessels = readKnossosRoi(out.root, out.prefix, bboxIn, 'uint32');
    % This line is to overwrite old results of detection(s) below
    vessels(vessels == 2) = 0; 
    vessels(vessels == 3) = 0;
    % Perform and myelin detection
    if  mean(bboxIn(3,:))>2000
        thresh= 70;
        else
        thresh= 60;
   end
    myelin = detectMyelinLocal(raw,thresh);
    %makeSegMovie(myelin, uint8(raw), '/home/alik/Desktop/test.avi');
    % Cut away additional (overlapping) context loaded to avoid border effects
    toRemove = abs(bboxIn - bboxOut);
    myelin = cutAwayBorder(myelin, toRemove);
    vessels = cutAwayBorder(vessels, toRemove);
    % Make sure old data is not overwritten
    myelin(vessels == 1) = 0;
    % Make sure no voxels are labeled double before pasting together
    vessels = uint32(vessels);
    vessels(myelin) = 3;
    gradientCorrected.root = '/gaba/u/alik/wKcubes/2017-04-03_ppcAK99_76x79_corrected/segmentation/1/';
    gradientCorrected.prefix = '2017-04-03_ppcAK99_76x79_corrected_mag1';
    writeKnossosRoi(gradientCorrected.root, gradientCorrected.prefix, bboxOut(:,1)', uint32(vessels), 'uint32', '', 'noRead');

end

function data = cutAwayBorder(data, border)

    data(1:border(1,1),:,:) = [];
    data(:,1:border(2,1),:) = [];
    data(:,:,1:border(3,1)) = [];
    data(end-border(1,2)+1:end,:,:) = [];
    data(:,end-border(2,2)+1:end,:) = [];
    data(:,:,end-border(3,2)+1:end) = [];

end

