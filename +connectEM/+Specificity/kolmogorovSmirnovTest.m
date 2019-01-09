function [pValue, ksStat] = ...
        kolmogorovSmirnovTest(testSamples, nullSamples, varargin)
    % [pValue, ksStat] = kolmogorovSmirnovTest(testSamples, nullSamples, varargin)
    %   Kolmogorov-Smirnov test that is slightly more versatile than
    %   MATLAB's built-in function. In particular, it allows weighted 
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.tail = 'equal';
    opt.nullWeights = ones(numel(nullSamples), 1);
    opt = Util.modifyStruct(opt, varargin{:});
    
   [nullKnots, ~, nullMasses] = unique(nullSamples(:));
    nullMasses = accumarray(nullMasses, opt.nullWeights);
    nullMasses = nullMasses / sum(nullMasses);
    
    file = tempname();
    write = @(name, val) writeFile(file, name, val);
    
    write('/null_knots', nullKnots);
    write('/null_masses', nullMasses);
    write('/observations', testSamples);
    
    % NOTE(amotta): We're shelling out to an R script for the Kolmogorov-
    % Smirnov test. The script is assume to have the same name as this
    % function with exception for the file suffix being `.r`.
    rScriptPath = strcat(mfilename('fullpath'), '.r');
    
    switch opt.tail
        case 'smaller', alternative = 'less';
        case 'equal', alternative = 'two.sided';
        case 'larger', alternative = 'greater';
        otherwise, error('Invalid tail option "%s"', opt.tail)
    end
    
    % HACK(amotta): Dammit MATLAB! MATLAB ships with its own GLIBC version
    % and accordingly modifies the library search path. As a result, all R
    % functions that were compiled against the system's GLIC will fail when
    % shelling out from MATLAB. Instead, we use `env -i` to start `Rscript`
    % in a clean environment without MATLAB's annoying search path.
    command = sprintf( ...
        'env -i Rscript %s %s %s', ...
        rScriptPath, alternative, file);
   [status, results] = system(command);
    
    assert(status == 0);
    results = strsplit(results);
    pValue = str2double(results{1});
    ksStat = str2double(results{2});
    
    % Clean up
    delete(file);
end

function writeFile(file, name, val)
    h5create(file, name, size(val));
    h5write(file, name, val);
end
