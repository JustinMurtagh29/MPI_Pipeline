function sampleDataset

    %load parameter
    load('~/results2/allParameter.mat');

    display('Sampling mean and standard deviation values for CNN normalization');
    
    % Specify how many 100^3 voxel cubes to sample
    nrCubesToSample = 200;
    sizeOfRoi = p.bbox(:,2) - p.bbox(:,1) + 1;
    
    MyVals = 0:255;
    accu = zeros(256,1);
    
    for i=1:nrCubesToSample
        lowerLeft = [randi(sizeOfRoi(1)-99); randi(sizeOfRoi(2)-99); randi(sizeOfRoi(3)-99)];
        lowerLeft = lowerLeft + p.bbox(:,1) - 1;
        bbox = cat(2,lowerLeft, lowerLeft + 99);
        raw = loadRawData(p.raw.root, p.raw.prefix, bbox, false);
        
        
        % Accumulate intensity values
        counts = histc(raw(:), MyVals);
        accu = accu + counts;
    end
    
    ClassicPDF = [accu'; MyVals];
    
    % Compute CDF
    N = sum(accu);
    cumulProb = zeros(1,256);
    tempCum = 0;
    for i = 1:256
        tempCum = tempCum + accu(i);
        cumulProb(1,i) = tempCum/N;
    end
    
    CDF = [cumulProb; MyVals];
    
    
    % Create a mapping function
    Mapper = fit(MyVals',cumulProb','linearinterp');
    save('~/SampleMy.mat', 'Mapper', 'ClassicPDF', 'CDF');
