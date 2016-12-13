classdef Identity < SynEM.Feature.TextureFeature
    %IDENTITY Wrapper for the identity feature.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    methods
        function obj = Identity()
            obj.name = 'Id';
            obj.border = [0, 0, 0];
            obj.numChannels = 1;
        end

        function raw = calc(~, raw)
            %Nothing to do
        end
    end
end
