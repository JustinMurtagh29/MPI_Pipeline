function asiT = apply(asiT)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % NOTE(amotta): This factor was measured with
    % connectEM.Consistency.Calibration.calibrateAsiAreas
    % git@gitlab.mpcdf.mpg.de:connectomics/pipeline.git c833deefadb65b63e8ababe56c07851d93bb398d (dirty)
    % amotta@m-01522. MATLAB 9.3.0.713579 (R2017b). 11-Feb-2019 15:22:35
    asiT.area = 0.825 * asiT.area;
end
