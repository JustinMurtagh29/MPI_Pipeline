function label = generateWeightLabels()

%% Feature names
label = cell(216,1);
label{1} = 'size border pixel';
label{2} = 'size smaller object';
label{3} = 'size bigger object';
label{4} = 'scalar products between directions';
label{5} = 'distance CoM';

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
                label{5+5*21*(i-1)+tempCounter*5+l} = [inputs{i} channels{j} 'size ' num2str(sizes{j}(k)) '.' num2str(k) ' ' stat{l}];
            end
            tempCounter = tempCounter + 1;
        end
    end
    tempCounter = 0;
end
label{224} = 'constant feature';



