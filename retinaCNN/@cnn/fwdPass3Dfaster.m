function fwdPass3Dfaster(cnet, input, result,bbox,normFunc)



%cnet.run.actvtTyp=@single;
% Load & normalize data and make sure single is used
activity = cell(cnet.numLayer, cnet.numFeature);
bboxWithBorder(:,1) = bbox(:,1) - ceil(cnet.randOfConvn'/2);
bboxWithBorder(:,2) = bbox(:,2) + ceil(cnet.randOfConvn'/2);
activity{1,1} = readKnossosRoi(input.root, input.prefix, bboxWithBorder);

if cnet.normalize
    activity{1,1} = normFunc(single(activity{1,1}));    
end

% Do the THING
layer = 2;
while size(activity,1) > 1
    for fm=1:cnet.layer{layer}.numFeature
        activity{2,fm} = zeros(size(activity{1,1}) - (cnet.filterSize - 1), class(single(1)));
        for oldFm=1:cnet.layer{layer-1}.numFeature
            activity{2, fm} = activity{2, fm} + convn(activity{1, oldFm}, ...
                single(cnet.layer{layer}.W{oldFm,fm}), 'valid');
        end
        activity{2, fm} = cnet.nonLinearity(activity{2, fm} + cnet.layer{layer}.B(fm));
    end
    activity(1,:) = [];
    layer = layer + 1;
end
% Save result to 3 KNOSSOS folders
activity(cellfun(@isempty,activity)) = [];
for i=1:size(activity,2)
    writeKnossosRoi([result.root result.folders{i} '/'], result.prefix, bbox(:,1)', single(activity{1,i}), 'single', '');
end
% Save result in single folder by taking the mean

writeKnossosRoi(result.root, result.prefix, bbox(:,1)', (single(activity{1,1})+single(activity{1,2})+single(activity{1,3}))/3, 'single', '');

end
