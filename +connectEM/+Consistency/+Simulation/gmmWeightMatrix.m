% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;
rng(0);

%% Configuration
preN = 100;
postN = 100;
synN = 20 * preN;
preClassN = 2;

mix = struct;
mix(1).name = 'Blue';
mix(1).mean = -0.484;
mix(1).std = 0.215;
mix(1).prob = 0.475;

mix(2).name = 'Red';
mix(2).mean = -0.909;
mix(2).std = 0.207;
mix(2).prob = 0.525;

outFile = '';

mix = struct2table(mix);
mix.prob = mix.prob / sum(mix.prob);
mix = sortrows(mix, 'prob', 'ascend');

info = Util.runInfo();
Util.showRunInfo(info);

%% Preparation
% We want the connection types to be distributed according to the
% probabilities in `mix.prob`.
assert(height(mix) == 2);
classProb = (1 + sqrt(2 * mix.prob(end) - 1)) / 2;

% Assign axons to axon classes
preClasses = linspace(0, 1, preN);
preClasses = discretize(preClasses, [0, classProb, 1]);
preClasses = preClasses(randperm(numel(preClasses)));

% Assign dendrites to target classes
postClasses = linspace(0, 1, postN);
postClasses = discretize(postClasses, [0, classProb, 1]);
postClasses = postClasses(randperm(numel(postClasses)));

flat = @(v) v(:);

%% Generate synapses
clear cur*;

asiT = table;
asiT.id = reshape(1:synN, [], 1);
asiT.preAggloId = flat(datasample(1:preN, synN, 'Replace', true));
asiT.postAggloId = flat(datasample(1:postN, synN, 'Replace', true));

asiT.connType = 1 + reshape( ...
    preClasses(asiT.preAggloId) ...
 == postClasses(asiT.postAggloId), [], 1);

asiT.log10Area = normrnd( ...
    mix.mean(asiT.connType), ...
    mix.std(asiT.connType));
asiT.area = 10 .^ asiT.log10Area;

% Add variables to compatibility with other scripts
asiT.type(:) = categorical({'PrimarySpine'});
asiT.axonClass(:) = categorical({'Corticocortical'});
asiT.targetClass(:) = categorical({'OtherDendrite'});

%% Save output
if ~isempty(outFile)
    out = struct;
    out.asiT = asiT;
    out.info = info;
    
    Util.saveStruct(outFile, out);
    Util.protect(outFile);
end
