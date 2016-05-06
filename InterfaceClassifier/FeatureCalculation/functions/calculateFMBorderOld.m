function border = calculateFMBorderOld( names, parameters, featureIndices )
%CALCULATEFMBORDER Calcualte the necessary border of the feature map.
% OUTPUT border: Maximal border of all filters in the feature map for each
%                dimension and for each direction. The total border for
%                each dimension is thus 2*border
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

border = zeros(1,3);
for feat = 1:length(parameters)
    if any(subsref(cell2mat(featureIndices{feat}),struct('type','()','subs',{{':'}})))
        parameter = parameters{feat};
        switch names{feat}
            case 'StructureTensor'
                filterBorder = parameter{2} + parameter{4};
                border = max(border,filterBorder);
            case {'GaussianSmoothed','DoG','LoG','GaussGrad','Hessian'}
                filterBorder = parameter{2};
                border = max(border,filterBorder);
            case {'Entropy','lPol','StdDev'}
                nhood = size(parameter{1});
                filterBorder = floor((nhood - 1)/2);
                border = max(border,filterBorder);
            case {'MaximumFilter','MinimumFilter'}
                nhood = parameter{1};
                if length(nhood) == 1
                    nhood = repmat(nhood,[1 3]);
                end
                filterBorder = floor((nhood - 1)/2);
                if length(parameter) == 2
                    sigma = parameter{2};
                    filterBorder = filterBorder + 2.*ceil(sigma*3);
                end
                border = max(border,filterBorder);
            case 'SphereAverage'
                r = parameter{1};
                filterBorder = [r,r,floor(r/2)];
                border = max(border,filterBorder);
        end
    end
end
border = border';

end
