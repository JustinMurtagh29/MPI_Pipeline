classdef Volume < SynEM.Feature.ShapeFeature
    %VOLUME Calculate the volume of a shape.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        func = [];
    end
    
    methods
        function obj = Volume(func, numFeat)
            %INPUT numFeat: (Optional) Number of features calculated by
            %           Volume.
            %           (Default: 1)
            %       func: Function handle that can be converted to string
            %           and that is applied to the volume of each shape.
            
            obj.name = 'Volume';
            if exist('numFeat','var') && ~isempty(numFeat)
                obj.numFeatures = numFeat;
            else
                obj.numFeatures = 1;
            end
            if exist('func','var') && ~isempty(func)
                if ischar(func)
                    obj.func = func;
                else
                    obj.func = func2str(func);
                end
            end
        end
        
        function feat = calculate(obj, vol, selFeat)
            % OUTPUT feat: [NxM] double with size = size(vol)
            
            if ~exist('selFeat','var') || isempty(selFeat)
                selFeat = true(obj.numFeatures,1);
            end
            vol = vol(:,selFeat);
            if isempty(obj.func)
                feat = cellfun(@(x)size(x,1),vol);
            else
                f = str2func(obj.func);
                feat = cellfun(@(x)f(size(x,1)),vol);
            end
        end
    end
    
end

