classdef EigenvaluesHessian < SynEM.Feature.TextureFeature
    %EIGENVALUESHESSIAN Eigenvalues of Hessian for 3d image data.
    % PROPERTIES
    % sigma: [3x1] array of double specifying the standard deviation
    %       of the gaussian kernel for each dimension. Each value must be bigger
    %       equal to zero and a zero corresponds to no filtering in the
    %       respective dimension.
    % truncate: Double. Filter is truncated at ceil(truncate.*sigma) many
    %       standard deviations.
    %       (Default: 3)
    % convMode: String specifying the convolution mode.
    %       (see convn).
    % sortMode: String specifying the sort mode. Options are
    %       'std': Sort eigenvalues in ascending order.
    %       'abs': Sort absolute values of eigenvalues in ascending order.
    %       (Legacy option for old feature calculation pipeline)
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        sigma
        truncate = 3;
        convMode = 'same';
        sortMode = 'std';
    end
    
    methods
        function obj = EigenvaluesHessian(sigma, truncate, convMode, ...
                sortMode)
            obj.name = 'EVsHessian';
            obj.sigma = sigma;
            if exist('truncate','var') && ~isempty(truncate)
                obj.truncate = truncate;
            end
            if exist('convMode','var') && ~isempty(convMode)
                obj.convMode = convMode;
            end
            if exist('sortMode','var') && ~isempty(sortMode)
                obj.sortMode = sortMode;
            end
            obj.numChannels = 3;
            obj.border = 2.*ceil(obj.sigma*obj.truncate);
        end
        
        function fm = calc(obj, raw)
            fm = obj.calculate(raw, obj.sigma, obj.truncate, ...
                obj.convMode, obj.sortMode);
        end
    end
    
    methods (Static)
        function fm = calculate(raw, sigma, truncate, convMode, sortMode)
            H = SynEM.Feature.Hessian.calculate(raw, sigma, truncate, ...
                convMode);
            [nx, ny, nz] = size(raw);
            ev = SynEM.Aux.eig3S([H{1,1}(:)';H{1,2}(:)';H{1,3}(:)'; ...
                H{2,2}(:)';H{2,3}(:)';H{3,3}(:)']);
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

