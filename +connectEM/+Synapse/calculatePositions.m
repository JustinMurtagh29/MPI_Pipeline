function pos = calculatePositions(param, syn, varargin)
    % NOTE(amotta): Factored out to auxiliary methods
    pos = Synapse.calculatePositions(param, syn.synapses, varargin{:});
end
