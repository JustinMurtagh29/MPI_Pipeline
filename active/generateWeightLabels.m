function label = generateWeightLabels()

%% Feature names
label = cell(224,1);
label{1} = 'surface between objects relative to size of the smaller object';
label{2} = 'relative size of objects';
label{3} = 'scalar product between major axes of objects';
label{4} = 'distance between objects';
label{5} = 'difference in bounding box edge size (dim 1)';
label{6} = 'difference in bounding box edge size (dim 2)';
label{7} = 'difference in bounding box edge size (dim 3)';
label{8} = 'difference in direction between objects (dim 1)';
label{9} = 'difference in direction between objects (dim 2)';
label{10} = 'difference in direction between objects (dim 3)';
label{11} = 'relative position of objects (dim 1)';
label{12} = 'relative position of objects (dim 2)';
label{13} = 'relative position of objects (dim 3)';

inputs = {'aff ' 'raw '};
channels = {'original ' 'sorted eigenvalues hessian ' 'gaussian smoothed gradient magnitude ' ...
    'intensity gaussian smoothed ' 'sorted eigenvalues structure ' 'laplace of gaussian ' 'difference of gaussians '};
sizes = {[1] [3 3 3 5 5 5] [3 5] [3 5] [3 3 3 5 5 5] [3 5] [3 5]};
stat = {'max' 'min' 'mean' 'median' 'standard deviation'};

tempCounter = 0;
for i=1:length(inputs)
    for j=1:length(channels)
        for k=1:length(sizes{j})
            for l=1:length(stat)
                label{13+5*21*(i-1)+tempCounter*5+l} = [inputs{i} channels{j} 'size ' num2str(sizes{j}(k)) '.' num2str(k) ' ' stat{l}];
            end
            tempCounter = tempCounter + 1;
        end
    end
    tempCounter = 0;
end
label{224} = 'constant feature';

end

