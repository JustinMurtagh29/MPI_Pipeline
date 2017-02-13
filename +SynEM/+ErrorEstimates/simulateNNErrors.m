function [ nnRP ] = simulateNNErrors( noNeurons, c_r, p, r, empDist, synT, runs )
%SIMULATENNERRORS Simulate synapse detection errors for connectome.
% INPUT noNeurons: Integer total number of neurons.
%       c_r: scalar double
%           connectivity ratio in connectome (must be < 1).
%       p: scalar double
%           Synapse detection precision.
%       r: scalar double
%           Synapse detection recall.
%       empDist: (Optional) [Nx1] double
%           Discrete probability distribution where the i-th entry
%           specified the probability that a connection has i synapses.
%           (Default: Feldmeyer dist)
%       synT: (Optional) scalar int
%           Number of synapses to accept a neuron-to-neuron connection as
%           connected.
%           (Default: 1)
%       runs: (Optional) scalar int
%           Number of simulation runs.
%           (Default: 1)
% OUTPUT nnRP: [runs x 2] double
%           Recall-precision pairs for the simulated synapse prediction.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('empDist','var') || isempty(empDist)
    empDist = [0 2 5 2 2]./11;
else %normalize to make sure
    empDist = empDist./sum(empDist);
end
if ~exist('synT','var') || isempty(synT)
    synT = 1;
end
if ~exist('runs','var') || isempty(runs)
    runs = 1;
end

nnRP = zeros(runs,2);
for i = 1:runs
    %create ground truth connectome
    cb = rand(noNeurons^2,1) < c_r; %binary connectome with ~c_r connections
    cw = zeros(length(cb),1);
    cw(cb) = discretesample(empDist, sum(cb(:))); %sample number of synapses

    %retrieve synapse with recall rate
    pred = zeros(length(cb),1);
    pred(cb) = binornd(cw(cb),r);

    %calculate number of fps
    numFPs = round((1-p)/p*r*sum(pred(:)));

    %distribute fps on connectome
    idx = randsample(length(pred),numFPs);
    fps = accumarray(idx,1);
    fps(end+1:length(pred)) = 0;
    pred = pred + fps;

    %calculate pr values
    tp = sum(cb & (pred >= synT));
    fp = sum(~cb & (pred >= synT));
    fn = sum(cb & (pred < synT));
    nnR = tp/(tp + fn);
    nnP = tp/(tp + fp);
    nnRP(i,:) = [nnR, nnP];
end
end

function x = discretesample(p, n)
% Samples from a discrete distribution
%
%   x = discretesample(p, n)
%       independently draws n samples (with replacement) from the 
%       distribution specified by p, where p is a probability array 
%       whose elements sum to 1.
%
%       Suppose the sample space comprises K distinct objects, then
%       p should be an array with K elements. In the output, x(i) = k
%       means that the k-th object is drawn at the i-th trial.
%       
%   Remarks
%   -------
%       - This function is mainly for efficient sampling in non-uniform 
%         distribution, which can be either parametric or non-parametric.         
%
%       - The function is implemented based on histc, which has been 
%         highly optimized by mathworks. The basic idea is to divide
%         the range [0, 1] into K bins, with the length of each bin 
%         proportional to the probability mass. And then, n values are
%         drawn from a uniform distribution in [0, 1], and the bins that
%         these values fall into are picked as results.
%
%       - This function can also be employed for continuous distribution
%         in 1D/2D dimensional space, where the distribution can be
%         effectively discretized.
%
%       - This function can also be useful for sampling from distributions
%         which can be considered as weighted sum of "modes". 
%         In this type of applications, you can first randomly choose 
%         a mode, and then sample from that mode. The process of choosing
%         a mode according to the weights can be accomplished with this
%         function.
%
%   Examples
%   --------
%       % sample from a uniform distribution for K objects.
%       p = ones(1, K) / K;
%       x = discretesample(p, n);
%
%       % sample from a non-uniform distribution given by user
%       x = discretesample([0.6 0.3 0.1], n);
%
%       % sample from a parametric discrete distribution with
%       % probability mass function given by f.
%       p = f(1:K);
%       x = discretesample(p, n);
%

%   Created by Dahua Lin, On Oct 27, 2008
%

%% parse and verify input arguments

assert(isfloat(p), 'discretesample:invalidarg', ...
    'p should be an array with floating-point value type.');

assert(isnumeric(n) && isscalar(n) && n >= 0 && n == fix(n), ...
    'discretesample:invalidarg', ...
    'n should be a nonnegative integer scalar.');

%% main

% process p if necessary

K = numel(p);
if ~isequal(size(p), [1, K])
    p = reshape(p, [1, K]);
end

% construct the bins

edges = [0, cumsum(p)];
s = edges(end);
if abs(s - 1) > eps
    edges = edges * (1 / s);
end

% draw bins

rv = rand(1, n);
c = histc(rv, edges);
ce = c(end);
c = c(1:end-1);
c(end) = c(end) + ce;

% extract samples

xv = find(c);

if numel(xv) == n  % each value is sampled at most once
    x = xv;
else	% some values are sampled more than once
    xc = c(xv);
    d = zeros(1, n);
    dv = [xv(1), diff(xv)];
    dp = [1, 1 + cumsum(xc(1:end-1))];
    d(dp) = dv;
    x = cumsum(d);
end

% randomly permute the sample's order
x = x(randperm(n));
end

