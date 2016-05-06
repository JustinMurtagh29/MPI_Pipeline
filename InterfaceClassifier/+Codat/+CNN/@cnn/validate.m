function loss = validate( cnet, stacksV, options )
%VALIDATE Calculate average pixel loss on input given the target.
% See cnet.train for documentation of inputs.
% NOTE If validation should be performed on GPU then the cnet must be
%      pushed to GPU before this function is called.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ischar(options.data_pre) %in case validate is called directly
    options.data_pre = str2func(options.data_pre);
end
loss = zeros(1,0,'single');
cnet = cnet.setConvMode(options.val_fwd_alg);
for s = 1:length(stacksV)
    [input, target, targetMask] = options.data_pre(stacksV{s});
    val_size = min(options.val_size,cnet.mSize(target,1:3)); %in case val_size is larger
    if options.gpuDev
        input = gpuArray(input);
        target = gpuArray(target);
    end
    iterations = prod(floor(cnet.mSize(target,1:3)./val_size));
    for iter = 1:iterations
        [x_test,y_test,curr_mask] = Codat.CNN.Misc.tileTrainingCubes(input, target, ...
                                    val_size, cnet.border, iterations, iter, targetMask);
        prediction = cnet.predict( x_test{1} );
        loss(end + 1) = gather(cnet.lossLayer( prediction, y_test{1}, curr_mask{1})); %#ok<AGROW>
    end
end
loss = sum(loss)/length(loss);

end
