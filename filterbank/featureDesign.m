function [weights weightNames] = featureDesign(input, border)

% Finding the pixeldata for each voxel of the borders
temp = arrayfun(@(b) real(input(b.PixelIdxList)), border, 'UniformOutput', false);

% Computing 7 features
% - 5 Quantiles: 0, 0.25, 0.5, 0.75, 1
weights = cell2mat(arrayfun(@(t) quantile(t{:}, [0, 0.25, 0.5, 0.75, 1])', temp, 'UniformOutput', false))';
% - Mean
weights = [weights arrayfun(@(t) mean(t{:}), temp)'];
% - Standard deviation
weights = [weights arrayfun(@(t) std(t{:}), temp)'];

weightNames = { 'quantile0', 'quantile0.25', 'quantile0.5', 'quantile0.75', 'quantile1', 'mean', 'std' };

end

