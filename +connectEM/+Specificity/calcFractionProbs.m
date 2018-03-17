function [fracs, probs] = calcFractionProbs(mnProbs, n)
    % [fracs, probs] = calcFractionProbs(mnProbs, n)
    %   Given the probabilities of three classes (A, B, and C) and a number
    %   of events (n), this function returns the probability distribution
    %   over the fraction A / (A + B) assuming a multinomial distribution
    %   with n trials.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % Sanity checks
    assert(isequal(size(mnProbs), [1, 3]));
    
    % Evaluate multinomial probability distribution.
    %
    % NOTE(amotta): The third entry becomes negative if the sum of the
    % first two entries becomes larger than `n`. This might seem bad. But
    % is in fact useful as it guarantees that `mnpdf` will return zero!
    grid = zeros(n + 1, n + 1, 3);
   [grid(:, :, 1), grid(:, :, 2)] = ndgrid(0:n, 0:n);
    grid(:, :, 3) = n - sum(grid, 3);

    grid = reshape(grid, [], 3);
    probs = mnpdf(grid, mnProbs);
    
    fracs = grid(:, 1) ./ sum(grid(:, 1:2), 2);
    fracs(isnan(fracs)) = 0;

   [fracs, ~, uniIds] = unique(fracs);
    probs = accumarray(uniIds, probs);
end