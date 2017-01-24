classdef EigenvaluesStructureTensor < SynEM.Feature.TextureFeature
    %EIGENVALUESSTRUCTURETENSOR Eigenvalues of structure tensor for 3d
    %image data.
    % PROPERTIES
    % sigmaD: Integer scalar or vector specifying the standard deviation
    %       of the gaussian derivative for all or seperately for each
    %       dimension of the input image.
    % sigmaW: Integer scalar or vector specifying the standard deviation
    %       of the gaussian window for all or seperately for each
    %       dimension of the input image.
    % truncateD: (Optional) Double. Truncate filter at this many standard
    %         deviations.
    %         (Default: 3)
    % truncateW: (Optional) Double. Truncate filter at this many standard
    %         deviations.
    %         (Default: 3)
    % convMode: String specifying the convolution mode.
    %       (see convn).
    % sortMode: String specifying the sort mode. Options are
    %       'std': Sort eigenvalues in ascending order.
    %       'abs': Sort absolute values of eigenvalues in ascending order.
    %       (Legacy option for old feature calculation pipeline)
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        sigmaD
        sigmaW
        truncateD = 3;
        truncateW = 3;
        convMode = 'same'
        sortMode = 'std';
    end
    
    methods
        function obj = EigenvaluesStructureTensor(sigmaD, sigmaW, ...
                truncateD, truncateW, convMode, sortMode)
            obj.name = 'EVsStructureTensor';
            obj.sigmaD = sigmaD;
            obj.sigmaW = sigmaW;
            if exist('truncateD','var') && ~isempty(truncateD)
                obj.truncateD = truncateD;
            end
            if exist('truncateW','var') && ~isempty(truncateW)
                obj.truncateW = truncateW;
            end
            if exist('convMode','var') && ~isempty(convMode)
                obj.convMode = convMode;
            end
            if exist('sortMode','var') && ~isempty(sortMode)
                obj.sortMode = sortMode;
            end
            obj.numChannels = 3;
            obj.border = 2.*(ceil(obj.sigmaW*obj.truncateW)  + ...
                ceil(obj.sigmaD*obj.truncateD));
        end
        
        function S = calc(obj, raw)
            S = obj.calculate(raw, obj.sigmaD, obj.sigmaW, ...
                obj.truncateD, obj.truncateW, obj.convMode, obj.sortMode);
        end
    end
    
    methods (Static)
        function fm = calculate(raw, sigmaD, sigmaW, ...
                truncateD, truncateW, convMode, sortMode)
            S = SynEM.Feature.StructureTensor.calculate(raw, sigmaD, ...
                sigmaW, truncateD, truncateW, convMode);
            [nx, ny, nz] = size(raw);
            ev = SynEM.Aux.eig3S([S{1,1}(:)';S{1,2}(:)';S{1,3}(:)'; ...
                        S{2,2}(:)';S{2,3}(:)';S{3,3}(:)']);
            switch sortMode
                case 'abs'
                    [~,sortIds] = sort(abs(ev),1);
                    linearIds = bsxfun(@plus, sortIds, ...
                        (0:(size(sortIds,2) - 1)).*3);
                    ev = ev(linearIds);
                    fm = cell(3,1);
                    fm{1} = reshape(ev(1,:),nx,ny,nz);
                    fm{2} = reshape(ev(2,:),nx,ny,nz);
                    fm{3} = reshape(ev(3,:),nx,ny,nz);
                case 'std'
                    fm{1} = reshape(ev(1,:),nx,ny,nz);
                    fm{2} = reshape(ev(2,:),nx,ny,nz);
                    fm{3} = reshape(ev(3,:),nx,ny,nz);
                otherwise
                    error('Unknown sort mode %s.',sortMode);
            end
        end
    end
    
end

