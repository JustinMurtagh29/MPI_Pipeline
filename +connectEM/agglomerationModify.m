function agglomerationModify(m, filename, graph, optional)

connectEM.agglomeration( ...
        m.borderSizeDendrites, m.segmentSizeDendrites, m.borderSizeAxons, m.segmentSizeAxons, ...
        m.axonProbThreshold, m.dendriteProbThreshold, m.spineProbThreshold, ...
        m.probThresholdDendrite, m.sizeThresholdDendrite, m.probThresholdAxon, m.sizeThresholdAxon, ...
        m.erProbThreshold, ...
        m.dendriteProbSpines, m.probThresholdSpines, m.maxStepsSpines, ...
        filename, graph, optional);
