classdef EVsHessian < SynEM.Feature.TextureFeature
    %EVSHESSIAN Eigenvalues for hessian sorted by increasing absolute
    % absolute value.
    % n_mean: The mean used for raw data normalization.
    % n_std: The standard deviation used for raw data normlaization.
    %
    % NOTE n_mean and n_std are used to ensure a consistent normalization
    %      w.r.t. to uint8 raw data.
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        sigma
        filterSiz
        n_mean = 0
        n_std = 1
    end

    methods
        function obj = EVsHessian(sigma, filterSiz, n_mean, n_std)
            obj.name = 'EVsHessian';
            obj.sigma = sigma;
            obj.numChannels = 3;
            obj.filterSiz = filterSiz;
            obj.border = 2.*filterSiz;
            if exist('n_mean', 'var') && ~isempty(n_mean)
                obj.n_mean = n_mean;
            end
            if exist('n_std', 'var') && ~isempty(n_std)
                obj.n_std = n_std;
            end
        end

        function fm = calc(obj, raw)
            
            %normalize
            raw = raw + obj.n_mean/obj.n_std;
            
            %calculate hessian
            sigma = obj.sigma; %#ok<*PROPLC>
            filterSiz = obj.filterSiz;
            sigma = num2cell(sigma);
            coords = arrayfun(@(siz)(-siz:siz).',filterSiz, ...
                'UniformOutput',false);
            coords = cellfun(@(coords,dim) ...
                permute(coords,[2:dim,1,(dim+1):3]), ...
                coords,num2cell(1:3),'UniformOutput',false);
            gauss = cellfun(@(coords,sigma)exp(-(coords./sigma).^2./2), ...
                coords,sigma,'UniformOutput',false);
            gauss = cellfun(@(gauss)gauss./sum(gauss),gauss, ...
                'UniformOutput',false);
            gaussD = cellfun(@(coords,gauss,sigma) ...
                -gauss.*coords./(sigma.^2),coords,gauss,sigma, ...
                'UniformOutput',false);
            gaussD2 = cellfun(@(coords,gauss,sigma) ...
                gauss.*(coords.^2./(sigma.^2)-1)./(sigma.^2), ...
                coords,gauss,sigma,'UniformOutput',false);
            H = cell(3,3);
            for dim1 = 1:3
                for dim2 = 1:dim1
                    if(dim1==dim2)
                        H{dim1,dim2} = convn(raw,gaussD2{dim1},'same');
                    else
                        H{dim1,dim2} = convn(raw,gaussD{dim1},'same');
                        H{dim1,dim2} = convn(H{dim1,dim2},gaussD{dim2}, ...
                                             'same');
                    end
                    for dim = setdiff(1:3,[dim1,dim2])
                        H{dim1,dim2} = convn(H{dim1,dim2},gauss{dim}, ...
                                             'same');
                    end
                end
            end
            [nx,ny,nz] = size(raw);
            ev = SynEM.Aux.eig3S([H{1,1}(:)'; H{2,1}(:)'; H{3,1}(:)'; ...
                                  H{2,2}(:)'; H{3,2}(:)'; H{3,3}(:)']);
            
            % sort by absolute value
            [~,sortIds] = sort(abs(ev),1);
            linearIds = bsxfun(@plus, sortIds, ...
                (0:(size(sortIds,2) - 1)).*3);
            ev = ev(linearIds);
            fm = cell(3,1);
            fm{1} = reshape(ev(1,:),nx,ny,nz);
            fm{2} = reshape(ev(2,:),nx,ny,nz);
            fm{3} = reshape(ev(3,:),nx,ny,nz);
        end
    end

end