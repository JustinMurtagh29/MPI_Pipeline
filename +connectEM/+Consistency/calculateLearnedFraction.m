function [learnedFrac, unlearnedFrac, cvThresh] = ...
        calculateLearnedFraction(synT, pairConfig, ctrlPairConfig, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.bandWidth = 0.06;
    opt.queryPoints = linspace(eps, sqrt(2), 101);
    opt.method = 'KernelDensity'; % for backward compatibility
    opt = Util.modifyStruct(opt, varargin{:});
    
    switch lower(opt.method)
        case 'kerneldensity', impl = @kernelDensityMethod;
        case 'maxdifference', impl = @maxDifferenceMethod;
        otherwise, error('Invalid method "%s"', opt.method);
    end
    
   [learnedFrac, unlearnedFrac, cvThresh] = ...
        impl(opt, synT, pairConfig, ctrlPairConfig);
end

function [learnedFrac, unlearnedFrac, cvThresh] = ...
        maxDifferenceMethod(~, synT, pairConfig, ctrlPairConfig)
    cv = @(vals) std(vals, 0, 2) ./ mean(vals, 2);
    pairCvs = sort(cv(synT.area(pairConfig.synIdPairs)), 'ascend');
    ctrlCvs = sort(cv(synT.area(ctrlPairConfig.synIdPairs)), 'ascend');
    
    t = table;
    t.cv = [pairCvs; ctrlCvs];
    t.pdf = [ ...
        +repelem(1 / numel(pairCvs), numel(pairCvs), 1); ...
        -repelem(1 / numel(ctrlCvs), numel(ctrlCvs), 1)];
    t = sortrows(t, 'cv', 'ascend');
    t.cdf = cumsum(t.pdf);
    
   [~, cvThreshIdx] = max(t.cdf);
    cvThresh = t.cv(cvThreshIdx);
    unlearnedFrac = mean(pairCvs > cvThresh);
    learnedFrac = t.cdf(cvThreshIdx);
end

function [learnedFrac, unlearnedFrac, cvThresh] = ...
        kernelDensityMethod(opt, synT, pairConfig, ctrlPairConfig)
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
    cvThresh = opt.queryPoints(threshId);
end

function cvDens = calculateCvDensity(opt, synT, pairConfig)
    cvs = synT.area(pairConfig.synIdPairs);
    cvs = std(cvs, 0, 2) ./ mean(cvs, 2);
    
    cvDens = ksdensity( ...
        cvs, opt.queryPoints, ...
        'Support', [-1E3, sqrt(2) + 1E3], ...
        'Bandwidth', opt.bandWidth, ...
        'BoundaryCorrection', 'reflection');
    cvDens = cvDens / sum(cvDens);
end
