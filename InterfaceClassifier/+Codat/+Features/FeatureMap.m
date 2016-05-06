classdef FeatureMap < handle
    %FEATUREMAP Voxel-based feature map (fm).
    %
    % PROPERTIES
    % features: Cell array of function handles to the features. Each feature
    %           function must accept the raw data as first parameter and an
    %           can have an arbitrary number of additional parameters.
    %           Use the method FeatureMap.addFeature to add new features to the
    %           feature map. The output of a feature function must be an array
    %           of the same size as the input array of a cell array where each
    %           cell contains an array of the same size of the input array
    %           (e.g. eigenvalues of structure tensor/hessian). Additional
    %           outputs are ignored.
    % parameters: Cell array where each cell contains the parameters for the
    %           respective feature handle in features as a cell array.
    % border: Integer array specifying the border for 'valid' convolutions of
    %         this feature map.
    % numFeatures: The total number of features per voxel in this feature map.
    %         This might differ from length(features) because some feature
    %         functions can potentially contribute multiple voxel features
    %         (e.g. eigenvalues of structure tensor/hessian).
    % counter: Counter to iterative over the feature via
    %          calculateNextFeature. The counter is automatically set to
    %          zero when the feature map is loaded from a mat-file.
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        features
        parameters
        border
        numFeatures
        counter
    end
    
    methods
        function obj = FeatureMap()
            obj.features = cell(0);
            obj.parameters = cell(0);
            obj.border = zeros(0,3);
            obj.numFeatures = 0;
            obj.counter = 0;
        end
        
        function addFeature(fm, feature, inputParams)
            %Add a feature to the feature map.
            % INPUT feature: Function handle for the feature. (If all
            %                features are currently saved as strings then this should
            %                also be passed as a string).
            %       inputParams: Cell array of input parameters fo feature.
            
            fm.features{end + 1} = feature;
            fm.parameters{end + 1} = inputParams;
            update(fm);
        end
        
        function update(fm)
            %Update border and numFeatures of feature map.
            borderTmp = [0, 0, 0];
            numFeaturesTmp = 0;
            for i = 1:length(fm.features)
                input = cat(2,{[]},fm.parameters{i});
                [~, bd, nF] = fm.features{i}(input{:});
                numFeaturesTmp = numFeaturesTmp + nF;
                if length(bd) == 1
                    bd = repmat(bd,[1, 3]);
                end
                borderTmp = max([borderTmp; bd],[],1);
            end
            fm.border = borderTmp;
            fm.numFeatures = numFeaturesTmp;
        end
        
        function deleteFeature(fm, featureIndices)
            %Delete all specified features.
            % INPUT featureIndices: Integer vector or logical array specifying
            %           the features to delete.
            fm.features(featureIndices) = [];
            fm.parameters(featureIndices) = [];
            update(fm);
        end
        
        function out = calculateFeature(fm, featIdx, raw)
            %Calculate the specified feature.
            % INPUT featIdx: Linear index or logical index for one feature.
            %       raw: The raw data
            out = fm.features{featIdx}(raw,fm.parameters{featIdx}{:});
        end
        
        function out = calculateNextFeature(fm, raw)
            %Increases the counter and calculates and calculates the next
            %feature.
            % INPUT raw: The raw data.
            
            fm.counter = fm.counter + 1;
            out = fm.calculateFeature(fm.counter, raw);
        end
        
        function continues = hasNext(fm)
            %Check if there are more features to iterate over.
            continues = fm.counter + 1 <= length(fm.features);
        end
        
        function setCounter(fm, c)
            %Set counter to specified value.
            if c == round(c) && c >= 0 && c <= length(fm.features)
                fm.counter = c;
            end
        end
        
        function sobj = saveobj(fm)
            % Save object in a struct. Handle to feature functions are
            % replaced by strings.
            
            sobj.features = cellfun(@func2str, fm.features, 'UniformOutput', false);
            sobj.parameters = fm.parameters;
            sobj.border = fm.border;
            sobj.numFeatures = fm.numFeatures;
        end
    end
    
    methods (Static)
        function obj = loadobj(s)
            % Load feature map from matfile.
            obj = Codat.Features.FeatureMap();
            obj.features = cellfun(@str2func, s.features, 'UniformOutput', false);
            obj.parameters = s.parameters;
            obj.border = s.border;
            obj.numFeatures = s.numFeatures;
        end
    end
    
end
