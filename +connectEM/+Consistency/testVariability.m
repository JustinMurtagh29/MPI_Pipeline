function [pVal, testStat] = testVariability(synT, pairConfig, ctrlConfig)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    pairCv = synT.area(pairConfig.synIdPairs);
    pairCv = std(pairCv, 0, 2) ./ mean(pairCv, 2);
    
    ctrlCv = synT.area(ctrlConfig.synIdPairs);
    ctrlCv = std(ctrlCv, 0, 2) ./ mean(ctrlCv, 2);
    
    if numel(pairCv) * numel(ctrlCv) ...
            / (numel(pairCv) + numel(ctrlCv)) < 4
        warning('Asymptotic p-value is inaccurate');
    end
    
   [testStat, pVal] = kstest2( ...
       pairCv, ctrlCv, 'tail', 'larger');
end
