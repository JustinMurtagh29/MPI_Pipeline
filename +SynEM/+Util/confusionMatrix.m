function [ C, precision, recall, fpr, specificity, accuracy, CI ] = confusionMatrix( classLabels, prediction )
%CONFUSIONMATRIX Confusion matrix for binary classification with additional
% performance measures.
% INPUT classLabels: [Nx1] Binary vector of true class labels or int array
%           specifying groups of inputs samples. For
%           each group it is enough to find one member to consider the
%           whole group found.
%       prediction: [Nx1] Binary vector of predictions
% OUTPUT C: Confusion matrix [tp, fn; fp, tn]
%        precision: Precision for class 1 (true)
%        recall: Recall (true positive rate) for class 1 (true)
%        fpr: Fall-out (false positive rate) for class 1 (true)
%        specificity: Specificity (true negative rate) for class 1 (true)
%        accuracy: Accuracy value.
%        CI: 2x2 cell array containing the indices of the input samples
%           for the confusion matrix entries. If a grouping is supplied
%           then only one representative for each detected/missed group is
%           stored in CI.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%confusion matrix entries
if islogical(classLabels)
    tp = sum(classLabels & prediction);
    fp = sum(~classLabels & prediction);
    fn = sum(classLabels & ~prediction);
    tn = sum(~classLabels & ~ prediction);
    
    if nargout > 6 % get indices for confusin matrix entries
        CI{1,1} = find(classLabels & prediction);
        CI{2,1} = find(~classLabels & prediction);
        CI{1,2} = find(classLabels & ~prediction);
        CI{2,2} = find(~classLabels & ~ prediction);
    end
elseif isnumeric(classLabels) && all(round(classLabels) == classLabels)
    group = classLabels;
    classLabels = group > 0;
    group = double(group);
    tp = length(setdiff(unique(group(prediction)),0));
    fp = sum(~classLabels & prediction);
    fn = length(setdiff(unique(group),0)) - tp;
    tn = sum(~classLabels & ~prediction);
    
    if nargout > 6 % get indices for confusion matrix entries
%         %get all interfaces for detected groups
%         tpI = setdiff(unique(group(prediction)),0);
%         CI{1,1} = find(ismember(group,tpI) & prediction);
        
        %get one representative for found groups
        tpI = prediction & group > 0;
        CI{1,1} = accumarray(Util.renumber(group(tpI)),find(tpI), ...
            [],@(x)x(1));
        
        %get one representative for each missed groups
        fnI = ismember(group,setdiff(unique(group(:)),[group(tpI);0]));
        CI{1,2} = accumarray(Util.renumber(group(fnI)),find(fnI), ...
            [],@(x)x(1));
                   
        CI{2,1} = find(~classLabels & prediction);
        CI{2,2} = find(~classLabels & ~ prediction);
    end
end

C = [tp,fn;fp,tn];

%performance measures
precision = tp/(tp + fp);
recall = tp/(tp + fn);
fpr = fp/(fp + tn);
specificity = tn/(fp + tn);
accuracy = (tp + tn)/(tp + fp + tn + fn);

end

