function [learnedFrac, unlearnedFrac, cvTresh] = ...
        calculateLearnedFraction(synT, pairConfig, ctrlPairConfig, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.bandWidth = 0.06;
    opt.queryPoints = linspace(eps, sqrt(2), 101);
    opt = Util.modifyStruct(opt, varargin{:});
    
    % Calculate densities
    pairDens = calculateCvDensity(opt, synT, pairConfig);
    ctrlDens = calculateCvDensity(opt, synT, ctrlPairConfig);
    
    % Determine optimal threshold
    threshId = sign(pairDens - ctrlDens);
   [~, threshId] = max(arrayfun( ...
        @(id) ...
            sum(threshId(1:id)) ...
          - sum(threshId((id + 1):end)), ...
        1:(numel(opt.queryPoints) - 1)));
    
    learnedFrac = sum(pairDens(1:threshId)) - sum(ctrlDens(1:threshId));
    unlearnedFrac = sum(pairDens((threshId + 1):end));
    cvTresh = opt.queryPoints(threshId);
end

function cvDens = calculateCvDensity(opt, synT, pairConfig)
    cvs = synT.area(pairConfig.synIdPairs);
    cvs = std(cvs, 0, 2) ./ mean(cvs, 2);
    
    cvDens = ksdensity( ...
        cvs, opt.queryPoints, ...
        'Support', [0, sqrt(2)], ...
        'Bandwidth', opt.bandWidth, ...
        'BoundaryCorrection', 'reflection');
    cvDens = cvDens / sum(cvDens);
end
