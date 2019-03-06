function skels = buildAggloAndFlightSkels(sagglos, baseSkel)
    % skels = buildAggloAndFlightSkels(sagglos, baseSkel)
    %   Splits super-agglomerates into the segment-based agglomerate, and
    %   flight parts, respectively. The components of a super-agglomerate
    %   are then collected in a skeleton object, each.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    if ~exist('baseSkel', 'var') || isempty(baseSkel)
        baseSkel = skeleton();
    end
    
   [aggloParts, flightParts] = ...
        Superagglos.splitIntoAgglosAndFlights(sagglos);
    
    skels = num2cell(size(sagglos));
    skels = repelem(baseSkel, skels{:});
    for curIdx = 1:numel(sagglos)
        curSkel = skels(curIdx);

        curAgglos = aggloParts{curIdx};
        curFlights = flightParts{curIdx};

        for curPartIdx = 1:numel(curAgglos)
            curSkel = curSkel.addTree( ...
                'Agglomerate', ...
                curAgglos(curPartIdx).nodes(:, 1:3), ...
                curAgglos(curPartIdx).edges, [0, 0, 1, 1]);
        end

        for curPartIdx = 1:numel(curFlights)
            curSkel = curSkel.addTree( ...
                'Flight', ...
                curFlights(curPartIdx).nodes(:, 1:3), ...
                curFlights(curPartIdx).edges, [1, 0, 0, 1]);
        end
        
        skels(curIdx) = curSkel;
    end
end
