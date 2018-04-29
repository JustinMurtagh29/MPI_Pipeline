function skel = buildNml(agglo, aggloId, skel)
    % skel = buildNml(agglo, skel)
    %   Builds an NML file which can be modified in webKNOSSOS and used to
    %   map back the manually made changes onto the super-agglomerate.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    if ~exist('skel', 'var') || isempty(skel)
        skel = skeleton();
    end
    
    treeName = sprintf('Agglomerate %d', aggloId);
    
    try
        description = skel.parameters.experiment.description;
        description = sprintf('%s. %s', description, treeName);
    catch
        % Looks like there was no description there yet.
        description = treeName;
    end
    
    comments = arrayfun( ...
        @num2str, 1:size(agglo.nodes, 1), ...
        'UniformOutput', false);
    comments = strcat({'Node '}, comments(:));
    
    skel = skel.addTree( ...
        treeName, agglo.nodes(:, 1:3), ...
        agglo.edges, [], [], comments);
    skel = skel.setDescription(description);
end
