function border = calculateFMBorder( names, parameters, featureIndices )
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
                [~,filterBorder] = Codat.Features.eigenvaluesOfStructureTensor([],parameter{:});
                border = max(border,filterBorder);
            case 'Hessian'
                [~,filterBorder] = Codat.Features.eigenvaluesOfHessian([],parameter{:});
                border = max(border,filterBorder);
            case 'GaussianSmoothed'
                [~,filterBorder] = Codat.Features.gaussianFilter([],parameter{:});
                border = max(border,filterBorder);
            case 'DoG'
                [~,filterBorder] = Codat.Features.differenceOfGaussians([],parameter{:});
                border = max(border,filterBorder);
            case 'LoG'
                [~,filterBorder] = Codat.Features.laplacianOfGaussian([],parameter{:});
                border = max(border,filterBorder);
            case 'GaussGrad'
                [~,filterBorder] = Codat.Features.gaussGradientMagnitude([],parameter{:});
                border = max(border,filterBorder);
            case 'SphereAverage'
                [~,filterBorder] = Codat.Features.localAverage([],parameter{:});
                border = max(border,filterBorder);
            case {'Entropy','MaximumFilter','MinimumFilter','lPol','StdDev'}
                nhood = parameter{1};
                if length(nhood) == 1
                    nhood = repmat(nhood,[1 3]);
                end
                filterBorder = floor((nhood - 1)/2);
                if length(parameter) == 2
                    sigma = parameter{2};
                    filterBorder = filterBorder + ceil(3.*sigma);
                end
                border = max(border,filterBorder);
        end
    end
end
border = border';

end
