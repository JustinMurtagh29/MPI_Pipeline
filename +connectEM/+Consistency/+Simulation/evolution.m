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
        'ltpLinearInf', 'ltpExponentialFinite', ...
        'ltdLinearZero', 'ltdLinearNonZero', ...
        'ltdExponentialZero', 'ltdExponentialNonZero'};
    assert(ismember({opts.method}, validMethods));
    
    method = str2func(opts.method);
    synAreas = method(opts, synAreas);
end

function synAreas = ltdLinearZero(opts, synAreas) %#ok
    mps = opts.mps;
    ns = opts.mins / mps;
    t = (0:ns) * mps;
    
    subRate = 0.5 / opts.mins;
    synAreas = synAreas - t(:) * subRate;
    synAreas(synAreas < 0) = nan;
end

function synAreas = ltdLinearNonZero(opts, synAreas) %#ok
    mps = opts.mps;
    ns = opts.mins / mps;
    t = linspace(0, 1, ns + 1);
    
    decayTarget = 10 ^ (-1.25);
    synAreas = (1 - t(:)) .* synAreas + t(:) * decayTarget;
end

function synAreas = ltdExponentialZero(opts, synAreas) %#ok
    mps = opts.mps;
    ns = opts.mins / mps;
    t = (0:ns) * mps;
    
    decayRate = log(0.5) * opts.mins / log(0.01);
    synAreas = synAreas .* (0.5 .^ (t(:) / decayRate));
end

function synAreas = ltdExponentialNonZero(opts, synAreas) %#ok
    mps = opts.mps;
    ns = opts.mins / mps;
    t = (0:ns) * mps;
    
    decayTarget = 10 ^ (-1.25);
    decayRate = log(0.5) * opts.mins / log(0.05);
    
    delta = synAreas - decayTarget;
    synAreas = synAreas - delta .* (1 - 0.5 .^ (t(:) / decayRate));
end

function synAreas = ltpLinearInf(opts, synAreas) %#ok
    mps = opts.mps;
    ns = opts.mins / mps;
    t = (0:ns) * mps;
    
    addRate = 0.5 / opts.mins;
    synAreas = synAreas + t(:) * addRate;
end

function synAreas = ltpExponentialFinite(opts, synAreas) %#ok
    mps = opts.mps;
    ns = opts.mins / mps;
    t = (0:ns) * mps;
    
    approachTarget = 10 ^ (-0.1);
    approachRate = log(0.5) * opts.mins / log(0.25);
    
    delta = approachTarget - synAreas;
    synAreas = synAreas + delta .* (1 - 0.5 .^ (t(:) / approachRate));
end
