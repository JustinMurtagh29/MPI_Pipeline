function synAreas = evolution(synAreas, varargin)
    % synAreas = evolution(synAreas, varargin)
    %   Simulates synapse area evolution under different LTP / LTD models.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.method = '';
    opts.mps = 0.25; % minutes per step
    opts.mins = 60; % minutes to simulate
    opts = Util.modifyStruct(opts, varargin{:});
    
    % It's a bit annoying that this list must be kept in synch with the
    % below functions. But it's the best compromise between convenience and
    % reliability that I see...
    validMethods = { ...
        'ltpAdditive', 'ltpSaturated', ...
        'ltdSubtractive', 'ltdDivisive', 'ltdSaturated'};
    assert(ismember({opts.method}, validMethods));
    
    method = str2func(opts.method);
    synAreas = method(opts, synAreas);
end

function synAreas = ltdSubtractive(opts, synAreas) %#ok
    mps = opts.mps;
    ns = opts.mins / mps;
    t = (0:ns) * mps;
    
    subRate = 1 / 30;
    synAreas = synAreas - t(:) * subRate;
    synAreas(synAreas < 0) = nan;
end

function synAreas = ltdDivisive(opts, synAreas) %#ok
    mps = opts.mps;
    ns = opts.mins / mps;
    t = (0:ns) * mps;
    
    decayRate = 5;
    synAreas = synAreas .* (0.5 .^ (t(:) / decayRate));
end

function synAreas = ltdSaturated(opts, synAreas) %#ok
    mps = opts.mps;
    ns = opts.mins / mps;
    t = (0:ns) * mps;
    
    decayRate = 5;
    decayTarget = 10 ^ (-1.25);
    
    delta = synAreas - decayTarget;
    synAreas = synAreas - delta .* (1 - 0.5 .^ (t(:) / decayRate));
end

function synAreas = ltpAdditive(opts, synAreas) %#ok
    mps = opts.mps;
    ns = opts.mins / mps;
    t = (0:ns) * mps;
    
    subRate = 1 / 30;
    synAreas = synAreas + t(:) * subRate;
end

function synAreas = ltpSaturated(opts, synAreas) %#ok
    mps = opts.mps;
    ns = opts.mins / mps;
    t = (0:ns) * mps;
    
    approachRate = 5;
    approachTarget = 10 ^ (-0.1);
    
    delta = approachTarget - synAreas;
    synAreas = synAreas + delta .* (1 - 0.5 .^ (t(:) / approachRate));
end
