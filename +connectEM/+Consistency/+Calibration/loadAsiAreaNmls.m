function calibT = loadAsiAreaNmls(param, nmlDir)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    nmlNames = dir(fullfile(nmlDir, '*.nml'));
    nmlNames([nmlNames.isdir]) = [];
    nmlNames = {nmlNames.name};
    
    calcArea = @connectEM.Consistency.Calibration.calculateAsiArea;
    
    calibT = [];
    for curId = 1:numel(nmlNames)
        curNmlName = nmlNames{curId};
        curNmlFile = fullfile(nmlDir, curNmlName);
        
        curParts = '^asi-(\d+)_axon-(\d+)_spine-head-(\d+)\.nml$';
        curParts = regexpi(curNmlName, curParts, 'tokens', 'once');
        curParts = cellfun(@str2double, curParts);
        
        curNml = slurpNml(curNmlFile);
        curArea = calcArea(param.raw.voxelSize, curNml);
        
        curCalib = struct;
        curCalib.asiId = curParts(1);
        curCalib.axonId = curParts(2);
        curCalib.shId = curParts(3);
        curCalib.area = curArea;
        curCalib.nmlFile = curNmlFile;
        
        % HACK(amotta): Lazy initialization
        if isempty(calibT); calibT = curCalib; end
        calibT(curId) = curCalib; %#ok
    end
    
    calibT = struct2table(calibT, 'AsArray', true);
end
