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
    
    edges = cat(1, testSamples(:), nullSamples(:));
    edges = cat(1, -inf, unique(edges), inf);
    
   [~, testCdf] = ismember(testSamples(:), edges);
   	testCdf = cumsum(accumarray(testCdf, 1, size(edges)));
    testCdf = testCdf / testCdf(end);
    
   [~, nullCdf] = ismember(nullSamples(:), edges);
    nullCdf = cumsum(accumarray(nullCdf, opt.nullWeights, size(edges)));
    nullCdf = nullCdf / nullCdf(end);
    
    switch opt.tail
        case 'smaller'
            ksStat = nullCdf - testCdf;
        case 'equal'
            ksStat = abs(testCdf - nullCdf);
        case 'larger'
            ksStat = testCdf - nullCdf;
        otherwise
            error('Invalid tail "%s"', opt.tail);
    end
    
    ksStat = max(ksStat(2:(end - 1)));
    if ~ksStat; pValue = 1; return; end
    
    pValue = ksStatToPValue( ...
        numel(testSamples), ksStat, ...
        ismember(opt.tail, {'smaller', 'larger'}));
end

%% Translate Kolmogorov-Smirnov statistic to p-value
% This part of extracted from MATLAB's kstest.
function pValue = ksStatToPValue(n, ksStat, tailed)
    % Copyright 1993-2013 The MathWorks, Inc.
    
    if ~(tailed)
        s = n*ksStat^2;

        % For d values that are in the far tail of the distribution (i.e.
        % p-values > .999), the following lines will speed up the computation
        % significantly, and provide accuracy up to 7 digits.
        if (s > 7.24) ||((s > 3.76) && (n > 99))
            pValue = 2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
        else
            % Express d as d = (k-h)/n, where k is a +ve integer and 0 < h < 1.
            k = ceil(ksStat*n);
            h = k - ksStat*n;
            m = 2*k-1;

            % Create the H matrix, which describes the CDF, as described in Marsaglia,
            % et al. 
            if m > 1
                c = 1./gamma((1:m)' + 1);

                r = zeros(1,m);
                r(1) = 1; 
                r(2) = 1;

                T = toeplitz(c,r);

                T(:,1) = T(:,1) - (h.^(1:m)')./gamma((1:m)' + 1);

                T(m,:) = fliplr(T(:,1)');
                T(m,1) = (1 - 2*h^m + max(0,2*h-1)^m)/gamma(m+1);
             else
                 T = (1 - 2*h^m + max(0,2*h-1)^m)/gamma(m+1);
             end

            % Scaling before raising the matrix to a power
            if ~isscalar(T)
                lmax = max(eig(T));
                T = (T./lmax)^n;
            else
                lmax = 1;
            end

            % Pr(Dn < d) = n!/n * tkk ,  where tkk is the kth element of Tn = T^n.
            % p-value = Pr(Dn > d) = 1-Pr(Dn < d)
            pValue = (1 - exp(gammaln(n+1) + n*log(lmax) - n*log(n)) * T(k,k));
        end
    else
        t = n * ksStat;
        k = ceil(t):n;       % k is the variable of summation
        pValue = sum(exp( log(t) - n*log(n)+ gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1) + ...
            k.*log(k - t) + (n-k-1).*log(t+n-k)));

    end
end
